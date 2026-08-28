# Handoff — SmartChemist × R-BRICS block-graph MCTS

Purpose: let a fresh session (Claude Code / VS Code) continue this project with
full context. Read this top-to-bottom once, then see **Next steps**.

---

## 1. Project goal

Build a **graph-based Monte Carlo Tree Search** that assembles drug-like
molecules from chemically meaningful blocks, guided by a **PUCT** prior (AlphaGo
style) whose `P(s,a)` comes from **conditional connection probabilities learned
from a real dataset** (MOSES, GuacaMol, or ZINC drug-like).

The blocks and their legal connections come from combining two algorithms:

- **SmartChemist** — hierarchical annotation of chemically/biologically
  meaningful substructure patterns (SMARTS). Smaller/simpler groups are
  *overshadowed* by larger/more relevant ones. Defines **what stays together**.
- **R-BRICS** (revised BRICS) — retrosynthetic bond breaking into synthons.
  Defines **where you are allowed to cut**.

The combination (the "hybrid joiner"): a molecule becomes a graph whose **nodes
are annotated motifs** and whose **edges are R-BRICS-breakable bonds *between*
motifs**. Over a corpus this yields `P(next_block | current_block, bond_tag)`,
which becomes the MCTS prior.

---

## 2. Repositories

| Repo | Role | Notes |
|------|------|-------|
| `EtienneReboul/smiles_blocks` | our project (MCTS + R-BRICS rewrite) | `smiles_blocks/` package; `mcts.py`, `retrosynthesis.py`, `rbrics_patterns.py` |
| `BiomedSciAI/r-BRICS` | upstream R-BRICS | reference only; `rBRICS_public.py` |
| `torbengutermuth/SmartChemist` | upstream annotator | Django app; SMARTS live in `smarts/*.csv` |

SmartChemist pattern CSVs (the annotation library, ~41k patterns):
`SmartChemist/smarts/functional_groups.csv` (~180),
`SmartChemist/smarts/cyclic.csv` (~40.7k whole-scaffold references),
`SmartChemist/smarts/biologicals.csv` (~65),
`SmartChemist/smarts/smarts_with_hierarchy.csv` (overshadow hierarchy by row index).

---

## 3. Delivered so far (two standalone modules, tested)

Both now live in `src/smart_mcts/` (this repo was scaffolded with the
[audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage)
template: src layout, `uv` + `justfile`, `pyproject.toml`, pytest, ruff/ty,
zensical docs, GitHub CI workflows). Set up with `uv sync`; run checks with
`just check` / `uv run pytest`.

Import as `from smart_mcts import smartchemist_annotator, hybrid_joiner`.
**Known gap:** `hybrid_joiner` still imports `smiles_blocks.rbrics_patterns` /
`smiles_blocks.retrosynthesis` from the separate `EtienneReboul/smiles_blocks`
repo — those modules were never delivered here, so the import fails until they
are vendored into `src/smart_mcts/` (or installed alongside).

### `smartchemist_annotator.py`
Django-free, fingerprint-screened re-implementation of SmartChemist matching.

- `SmartChemistAnnotator.from_csvs(functional_groups=, cyclic=, biologicals=, hierarchy_csv=)`
  compiles the library once (~30 s for 41k patterns). `.save(path)` / `.load(path)`
  cache it (pickle) so you pay the compile cost once.
- `.annotate(mol, remove_overshadowed=False)` → list of match dicts:
  `{atom_indices, name, group, smarts, n_bonds, hierarchy, index, overshadowed}`.
  Output shape mirrors the original `SmartChemist._match_smarts_patterns`.
- Per molecule: cheap descriptor prefilter → RDKit **pattern-fingerprint subset
  screen** (`DataStructs.AllProbeBitsMatch`) → `HasSubstructMatch` on survivors.
- **Measured**: identical results to the naive full-library scan; **5.1×**
  faster (179 → 35 ms/mol over the full 41k library); screen leaves ~0.25%
  survivors. Speedup grows with library size (survivors stay ~constant).
- **Bug fixed vs upstream**: upstream `_check_overshadowed_patterns` had
  `if bonds1 > bonds2 ... elif bonds2 < bonds1` — the `elif` repeats the same
  condition and is dead code. Should be `bonds2 > bonds1`. Fixed here in
  `_resolve_overshadowing`.

### `hybrid_joiner.py`  ← this is the novel core
Combines annotator output + R-BRICS bonds into the block graph.

- `find_rbrics_bonds_live(mol, rbrics=None)` → `[RBond(a, b, tag_a, tag_b)]` in
  **live** atom indices (reuses `smiles_blocks.rbrics_patterns` dataclasses).
  NB: this differs from the repo's `retrosynthesis.annotate_rbond`, which keys
  bonds by **canonical rank** — those don't align with SmartChemist's live atom
  indices, so the joiner uses live indices throughout.
- `build_block_graph(mol, annotations, rbrics_bonds=None, *, tag_registry,
  unannotated_policy="absorb", motif_integrity=True)` → `BlockGraph`.
- `BlockGraph.blocks` = `Block(block_id, atom_indices, motif_names, smiles, region_ids)`;
  `.edges` = `BlockEdge(block_i, block_j, atom_i, atom_j, tag_i, tag_j)`.
- `BlockGraph.transition_records()` → directed, tag-annotated rows
  (both directions per edge) for corpus counting.
- `reassemble_canonical(mol, bg)` → rebuilds the molecule from blocks+edges;
  **verified lossless** (== original canonical SMILES) on the whole test panel.

**Block SMILES carry dummy atoms whose isotope = the *partner's* R-BRICS tag id**
(BRICS convention: `[16*]` = "was attached to an L16 environment"), via a
corpus-wide `TagRegistry`. This makes blocks self-describing for reassembly and
compatibility lookups. Verified `FragmentOnBonds` dummy-label semantics:
`dummyLabels=[(beginLabel, endLabel)]` puts `endLabel` on the begin-atom's
fragment and `beginLabel` on the end-atom's fragment.

---

## 4. Key design decisions (already made, revisit if needed)

**Hybrid rule.** Blocks = connected components after deleting *only* R-BRICS
bonds whose endpoints are in different motif regions. Consequences: an R-BRICS
bond inside one motif is never cut (the motif survives whole); two motifs joined
by a non-R-BRICS bond stay merged.

**`motif_integrity=True` (default).** A winning motif is never split by an
internal R-BRICS bond; overlapping motifs *merge* into one block (union-find over
motif atoms). This is the faithful reading of "annotated patterns are the nodes."
Discovered necessity from testing: without it, the amide C–N bond was cut in
paracetamol (an overlapping Acetyl motif stole the carbonyl) but kept whole in
biphenyl-amide — the same bond treated inconsistently. `False` = more R-BRICS-like
(overlapping motifs may fragment each other), denser vocabulary.

**`unannotated_policy` for linker atoms** (three genuinely distinct behaviors,
all lossless):
- `absorb` (default): linkers join the nearest motif region (multi-source BFS).
- `linker`: each connected unannotated region becomes its own block (finest
  granularity, explicit linker blocks).
- `merge`: only ever cut bonds between two *annotated* motifs (linkers never
  separated).

---

## 5. Findings on the existing repo code

**`retrosynthesis.py` / `rbrics_patterns.py` (the R-BRICS rewrite): good.**
Cleaner than upstream `FindrBRICSBonds` (precompiled patterns, `check_chemical_group`
prefilter, `RBRICSCompatibilityMap` replaces opaque `reactionDefs`). **Caveat:**
it keys bonds by `CanonicalRankAtoms`, which is **not stable across edits** to an
`RWMol`. Fine for storage/dedup, but keep live atom indices alongside if you ever
map a stored bond back onto a growing molecule. The joiner sidesteps this by
working in live indices.

**`smart_chemist.py`: the bottleneck.** Loops `AnnotatedPattern.objects.all()`
(41k rows) with a `HasSubstructMatch` per row **per molecule** — fatal at dataset
scale. Also `_check_overshadowed_patterns` is O(n²) and has the dead-`elif` bug
above. Superseded by `smartchemist_annotator.py`.

**`mcts.py`: already heavily optimized** (precomputed CDFs, `__slots__`, O(1)
`clone()`, virtual loss, shared-memory Arrow transport). **Two gaps for this
project:**
1. **Assembly is string concatenation** — `assembled_smiles` does
   `"".join(rows_by_uid[uid]["block"] for uid in chain)`. That only works for
   concatenable SMILES substrings, not R-BRICS blocks with dummy atoms, which
   must be reconnected at matching attachment points (`Chem.molzip` / BRICS
   reverse reactions on an `RWMol`).
2. **`MolState` is a *linear chain*** (`_uid_chain`, single `current_end_tag`,
   one-in-one-out). The joiner produces **branched** blocks with 1–4 attachment
   points (e.g. a benzene core with 3 substituents; imatinib's central rings).
   A linear chain cannot represent these. **This is the biggest architectural
   gap for the stated "graph-based MCTS" goal:** either restrict the vocabulary
   to ≤2 attachment points, or extend `MolState` to a real graph where the state
   tracks a *frontier of open attachment points* and PUCT ranges over
   `(open_point, compatible_block)` pairs. The conditional table would then key
   on `(block, open_point_tag)`.

---

## 6. Conditional-probability → PUCT mapping

`mcts.py` already implements the AlphaGo mapping: `FragmentLibrary._build_conditional_cache`
turns corpus `proba` into the `P(s,a)` prior used in `MCTSNode.puct`. To feed it
from the hybrid corpus:

- **Transition key should be `(block_a_id, bond_tag) → block_b_id`**, not just
  `block_a → block_b`. The same two motifs can join via different R-BRICS
  chemistry; conditioning on the tag keeps the prior chemically honest.
- **`unique_id` = canonical block SMILES-with-dummies** (optionally + motif
  label). `begin_tag`/`end_tag` come from the dummy tags.
- **Sparsity / backoff (do this from the start):** with ~41k motifs the specific
  counts are mostly 0/1. Emit three levels so you can interpolate (Katz /
  Kneser-Ney style): `P(b | a, tag)` → `P(b | motif_class(a), tag)` → `P(b | tag)`.

FragmentLibrary expected schemas (from `mcts.py`):
- fragment table columns: `unique_id, block, begin_tag, end_tag, frequency`
- conditional table columns: `unique_id, next_unique_id, frequency, proba`

---

## 7. Next steps (priority order)

1. **Corpus counting pipeline (thread 2).** Run annotator → joiner over a MOSES
   / GuacaMol subset; aggregate `BlockGraph.transition_records()` into the
   `(unique_id, next_unique_id, frequency, proba)` table **with the three backoff
   levels**. Decide `motif_integrity` and `unannotated_policy` *before* the full
   run (they change the vocabulary). Parallelize with multiprocessing; build the
   annotator index once and `.load()` per worker.
2. **Scale-test the joiner first** on a few thousand real MOSES molecules to
   measure block-vocabulary size and catch fragmentation edge cases before
   committing to a full corpus run. Cheap insurance.
3. **Graph assembly for the MCTS (thread 3 + the architectural gap in §5).**
   Replace string-concat with `RWMol`/`molzip` assembly, and extend `MolState`
   from a linear chain to a frontier-of-attachment-points graph so branched
   blocks are representable. These two are coupled; do them together.
4. **Fold the annotator speedups further** if corpus runtime demands it: replace
   the per-pattern `AllProbeBitsMatch` loop with a single vectorized numpy
   bit-matrix AND; add an exact canonical-SMILES hash fast-path for the 40k
   whole-scaffold "cyclic" references; process-parallelize.
5. **Upstream courtesy:** the SmartChemist overshadow bond bug (§3) is worth a PR.

---

## 8. Reproduce the validation

```python
from rdkit import Chem
from smartchemist_annotator import SmartChemistAnnotator
from hybrid_joiner import (build_block_graph, find_rbrics_bonds_live,
                           prepare_rbrics, reassemble_canonical, TagRegistry)

# build + cache the annotator index once
ann = SmartChemistAnnotator.from_csvs(
    functional_groups="SmartChemist/smarts/functional_groups.csv",
    cyclic="SmartChemist/smarts/cyclic.csv",
    biologicals="SmartChemist/smarts/biologicals.csv",
    hierarchy_csv="SmartChemist/smarts/smarts_with_hierarchy.csv",
)
ann.save("smartchemist_index.pkl")

rbrics, reg = prepare_rbrics(), TagRegistry()
mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")           # paracetamol
bg = build_block_graph(mol, ann.annotate(mol),
                       find_rbrics_bonds_live(mol, rbrics),
                       tag_registry=reg, motif_integrity=True)
assert reassemble_canonical(mol, bg) == Chem.MolToSmiles(mol)   # lossless
for b in bg.blocks: print(b.block_id, b.smiles, b.motif_names)
for e in bg.edges:  print(e.block_i, e.block_j, e.tag_i, e.tag_j)
```

Expected paracetamol (motif_integrity=True): 2 blocks —
`[2*]NC(C)=O` (Acetyl+Amide, whole) and `[4*]c1ccc(O)cc1` (Phenol) — one
`L5~L16` edge.

---

## 9. Open questions for the user

- `motif_integrity` on (faithful, sparser vocab) or off (R-BRICS-like, denser)?
- Which `unannotated_policy` — `absorb` (compact) or `linker` (explicit linker
  blocks, often better for generation)?
- Which dataset first, and full ~41k motif vocabulary or a curated subset?
- Are you committing to branched graph-MCTS (§5 gap), or restricting the block
  vocabulary to linear/≤2-attachment blocks for a first pass?
