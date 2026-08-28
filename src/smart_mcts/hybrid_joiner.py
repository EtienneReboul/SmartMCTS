"""Hybrid SmartChemist x R-BRICS fragmentation joiner.

Combines two decompositions of a molecule:

  * SmartChemist annotations  -> *what stays together* (chemically/biologically
    meaningful motifs, as atom-index sets).
  * R-BRICS breakable bonds   -> *where you are allowed to cut* (retrosynthetic
    joints, as atom pairs carrying environment tags like L16, L5, ...).

Hybrid rule
-----------
A block is a connected component of the molecular graph after deleting *only*
the R-BRICS bonds whose two endpoints fall in different motif regions. Three
consequences fall out for free:

  * an R-BRICS bond *inside* a single motif (e.g. an amide inside an annotated
    peptide) is never cut -> the motif survives whole;
  * two motifs joined by a bond R-BRICS would not break stay in one block;
  * unannotated linker atoms are handled by a configurable policy (default:
    absorbed into the nearest motif by graph distance, which splits a
    linker between the two motifs it bridges).

The output `BlockGraph` is self-describing and reassemblable: every block is a
sub-SMILES whose dummy atoms carry the *partner's* R-BRICS tag id (BRICS
convention), and every edge records the (tag_i, tag_j) pair driving the
conditional-probability model.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

from rdkit import Chem

# Reuse the user's curated R-BRICS pattern definitions.
from smart_mcts.rbrics_patterns import (
    RBRICSChemicalGroups,
    RBRICSRetrosynthesisBound,
)


# ===================================================================== #
#  Tag registry  (stable string-tag <-> int isotope label, corpus-wide) #
# ===================================================================== #
class TagRegistry:
    """Maps R-BRICS tags ('L16', 'L14b', ...) to stable small int ids.

    Dummy-atom isotopes must be ints; the *same* registry must be used across
    the whole corpus so a block's `[id*]` attachment point means the same thing
    in every molecule.
    """

    def __init__(self) -> None:
        self._t2i: dict[str, int] = {"no_tag": 0}
        self._i2t: dict[int, str] = {0: "no_tag"}

    def id(self, tag: str) -> int:
        if tag not in self._t2i:
            i = len(self._t2i)
            self._t2i[tag] = i
            self._i2t[i] = tag
        return self._t2i[tag]

    def tag(self, i: int) -> str:
        return self._i2t.get(i, "no_tag")

    def as_dict(self) -> dict[str, int]:
        return dict(self._t2i)


# ===================================================================== #
#  R-BRICS bond detection in LIVE atom indices                          #
# ===================================================================== #
@dataclass
class RBond:
    """One R-BRICS-breakable bond in live atom indices."""

    a: int
    b: int
    tag_a: str
    tag_b: str


def prepare_rbrics() -> dict:
    """Compile R-BRICS chemical-group and bond patterns once."""
    groups = {
        tag: Chem.MolFromSmarts(p) for tag, p in RBRICSChemicalGroups().patterns.items()
    }
    bonds = {}
    for tag, p in RBRICSRetrosynthesisBound().patterns.items():
        t1, t2 = tag.split("-")
        q = Chem.MolFromSmarts(p)
        if q is not None:
            bonds[tag] = (t1, t2, q)
    return {"groups": groups, "bonds": bonds}


def find_rbrics_bonds_live(mol: Chem.Mol, rbrics: dict | None = None) -> list[RBond]:
    """Return R-BRICS-breakable bonds with **live** atom indices and tags.

    Ring bonds are skipped (cannot be safely cleaved in a SMILES block).
    Deduplicated on the unordered atom pair. Unlike the canonical-rank version,
    these indices align directly with SmartChemist atom_indices.
    """
    rbrics = rbrics or prepare_rbrics()
    groups, bonds = rbrics["groups"], rbrics["bonds"]

    # cheap prefilter: which environments exist at all
    present = {t for t, q in groups.items() if q is not None and mol.HasSubstructMatch(q)}

    seen: set[frozenset] = set()
    out: list[RBond] = []
    for _label, (t1, t2, q) in bonds.items():
        if t1 not in present or t2 not in present:
            continue
        for a, b in mol.GetSubstructMatches(q):
            key = frozenset((a, b))
            if key in seen:
                continue
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond is None or bond.IsInRing():
                continue
            seen.add(key)
            out.append(RBond(a=a, b=b, tag_a=t1, tag_b=t2))
    return out


# ===================================================================== #
#  Block graph data structures                                          #
# ===================================================================== #
@dataclass
class Block:
    block_id: int
    atom_indices: tuple[int, ...]        # original (live) heavy-atom indices
    motif_names: tuple[str, ...]         # distinct winning motifs covering it
    smiles: str                          # sub-SMILES with [id*] dummy attachments
    region_ids: tuple[int, ...] = ()     # internal partition regions merged here


@dataclass
class BlockEdge:
    block_i: int
    block_j: int
    atom_i: int                          # original index in block_i
    atom_j: int                          # original index in block_j
    tag_i: str                           # R-BRICS env tag on the block_i side
    tag_j: str                           # R-BRICS env tag on the block_j side


@dataclass
class BlockGraph:
    smiles: str
    blocks: list[Block] = field(default_factory=list)
    edges: list[BlockEdge] = field(default_factory=list)

    def transition_records(self) -> list[dict]:
        """Directed (block -> block, tagged) records for corpus counting.

        Each undirected edge yields both directions, so frequencies are not
        artificially halved when you aggregate P(next | current, tag).
        """
        recs = []
        by_id = {b.block_id: b for b in self.blocks}
        for e in self.edges:
            bi, bj = by_id[e.block_i], by_id[e.block_j]
            recs.append(
                {
                    "from_smiles": bi.smiles,
                    "to_smiles": bj.smiles,
                    "from_motif": bi.motif_names,
                    "to_motif": bj.motif_names,
                    "tag_from": e.tag_i,
                    "tag_to": e.tag_j,
                }
            )
            recs.append(
                {
                    "from_smiles": bj.smiles,
                    "to_smiles": bi.smiles,
                    "from_motif": bj.motif_names,
                    "to_motif": bi.motif_names,
                    "tag_from": e.tag_j,
                    "tag_to": e.tag_i,
                }
            )
        return recs


# default specificity ranking for resolving overlapping motifs
_GROUP_RANK = {"biological": 0, "cyclic": 1, "functional_group": 2}


def _assign_regions(
    mol: Chem.Mol,
    winning: list[dict],
    policy: str,
    group_rank: dict[str, int],
    motif_integrity: bool,
) -> list[int]:
    """Assign every atom an integer region id.

    If `motif_integrity` (default), every winning motif is forced whole: all of
    its atoms are union-found into one region, so an internal R-BRICS bond is
    never cut and two *overlapping* motifs merge into a single block. This is
    the faithful reading of "annotated patterns connected by breakable bonds".

    If not, motif atoms take their highest-priority motif's region (overlapping
    motifs may fragment each other at R-BRICS joints — more R-BRICS-like).

    Unannotated atoms are then handled by `policy`:
        'absorb' : multi-source BFS -> nearest motif region (ties: lower rank)
        'linker' : each connected unannotated region gets its own fresh id
        'merge'  : left unassigned (-1) so bonds touching them are never cut
    """
    n = mol.GetNumAtoms()
    region = [-1] * n

    if motif_integrity:
        # union-find over atoms; union all atoms within each winning motif
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)

        for w in winning:
            atoms = w["atom_indices"]
            for at in atoms[1:]:
                union(atoms[0], at)
        # relabel roots of motif atoms to compact region ids
        roots: dict[int, int] = {}
        for w in winning:
            for at in w["atom_indices"]:
                r = find(at)
                if r not in roots:
                    roots[r] = len(roots)
                region[at] = roots[r]
    else:
        order = sorted(
            range(len(winning)),
            key=lambda k: (
                group_rank.get(winning[k]["group"], 99),
                -len(winning[k]["atom_indices"]),
                -winning[k].get("n_bonds", 0),
                winning[k].get("index", 0),
            ),
        )
        motif_region = {mk: pos for pos, mk in enumerate(order)}
        for mk in order:                      # assign low priority first ...
            for at in winning[mk]["atom_indices"]:
                region[at] = motif_region[mk]
        for mk in reversed(order):            # ... then overwrite with high priority
            for at in winning[mk]["atom_indices"]:
                region[at] = motif_region[mk]

    if policy == "merge":
        return region

    if policy == "absorb":
        # multi-source BFS from all motif atoms simultaneously
        dq = deque(i for i in range(n) if region[i] >= 0)
        while dq:
            i = dq.popleft()
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if region[j] == -1:
                    region[j] = region[i]
                    dq.append(j)

    # remaining -1 atoms (unreachable from any motif, or policy == 'linker'):
    # give each connected unannotated component a fresh region id
    next_id = (max(region) + 1) if region and max(region) >= 0 else 0
    for i in range(n):
        if region[i] != -1:
            continue
        comp = []
        stack = [i]
        region[i] = next_id
        while stack:
            x = stack.pop()
            comp.append(x)
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                j = nb.GetIdx()
                if region[j] == -1:
                    region[j] = next_id
                    stack.append(j)
        next_id += 1
    return region


def build_block_graph(
    mol: Chem.Mol,
    annotations: list[dict],
    rbrics_bonds: list[RBond] | None = None,
    *,
    tag_registry: TagRegistry | None = None,
    unannotated_policy: str = "absorb",
    motif_integrity: bool = True,
    group_rank: dict[str, int] | None = None,
    rbrics: dict | None = None,
) -> BlockGraph:
    """Build the hybrid block graph for one molecule.

    Parameters
    ----------
    mol : an RDKit mol (no explicit Hs needed).
    annotations : SmartChemist match dicts (each with 'atom_indices', 'group',
        'name', optional 'n_bonds'/'index'/'overshadowed'). Overshadowed
        matches are ignored.
    rbrics_bonds : precomputed live-index R-BRICS bonds; computed if None.
    tag_registry : shared TagRegistry (build one per corpus run, not per mol).
    unannotated_policy : 'absorb' | 'linker' | 'merge'.
    motif_integrity : if True (default), winning motifs are never split and
        overlapping motifs merge into one block.
    """
    reg = tag_registry or TagRegistry()
    group_rank = group_rank or _GROUP_RANK
    if rbrics_bonds is None:
        rbrics_bonds = find_rbrics_bonds_live(mol, rbrics)

    winning = [a for a in annotations if not a.get("overshadowed", False)]
    region = _assign_regions(
        mol, winning, unannotated_policy, group_rank, motif_integrity
    )

    # which R-BRICS bonds straddle two regions -> these are the cuts.
    # Under 'merge', a bond touching an unannotated atom (region -1) is never a
    # cut, so only genuinely inter-motif bonds separate blocks.
    cut_specs: list[tuple[int, RBond]] = []
    for rb in rbrics_bonds:
        ra, rb_ = region[rb.a], region[rb.b]
        if ra == rb_:
            continue
        if unannotated_policy == "merge" and (ra < 0 or rb_ < 0):
            continue
        bond = mol.GetBondBetweenAtoms(rb.a, rb.b)
        if bond is not None:
            cut_specs.append((bond.GetIdx(), rb))

    smiles = Chem.MolToSmiles(mol)

    # No cuts -> whole molecule is a single block
    if not cut_specs:
        motifs = tuple(sorted({w["name"] for w in winning}))
        return BlockGraph(
            smiles=smiles,
            blocks=[Block(0, tuple(range(mol.GetNumAtoms())), motifs, smiles)],
            edges=[],
        )

    bond_idxs = [bi for bi, _ in cut_specs]
    # dummyLabels=(beginLabel, endLabel): the begin-side fragment receives
    # endLabel, the end-side fragment receives beginLabel -> each block's dummy
    # carries the PARTNER tag. Map each cut bond's (a,b) to begin/end correctly.
    dummy_labels = []
    for bi, rb in cut_specs:
        bond = mol.GetBondWithIdx(bi)
        if bond.GetBeginAtomIdx() == rb.a:
            dummy_labels.append((reg.id(rb.tag_a), reg.id(rb.tag_b)))
        else:
            dummy_labels.append((reg.id(rb.tag_b), reg.id(rb.tag_a)))

    frag_mol = Chem.FragmentOnBonds(
        mol, bond_idxs, addDummies=True, dummyLabels=dummy_labels
    )
    frag_atom_tuples = Chem.GetMolFrags(frag_mol, asMols=False, sanitizeFrags=False)
    frag_mols = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)

    # map original atom idx -> block id (dummies have idx >= original N)
    n_orig = mol.GetNumAtoms()
    atom2block: dict[int, int] = {}
    blocks: list[Block] = []
    motif_by_atom: dict[int, set[str]] = {}
    for w in winning:
        for at in w["atom_indices"]:
            motif_by_atom.setdefault(at, set()).add(w["name"])

    for blk_id, (atom_tuple, fmol) in enumerate(zip(frag_atom_tuples, frag_mols)):
        orig_atoms = tuple(a for a in atom_tuple if a < n_orig)
        for a in orig_atoms:
            atom2block[a] = blk_id
        motifs = sorted({m for a in orig_atoms for m in motif_by_atom.get(a, ())})
        regions = tuple(sorted({region[a] for a in orig_atoms}))
        try:
            blk_smiles = Chem.MolToSmiles(fmol)
        except Exception:
            blk_smiles = Chem.MolToSmiles(fmol, canonical=False)
        blocks.append(
            Block(blk_id, orig_atoms, tuple(motifs), blk_smiles, regions)
        )

    edges = []
    for _bi, rb in cut_specs:
        edges.append(
            BlockEdge(
                block_i=atom2block[rb.a],
                block_j=atom2block[rb.b],
                atom_i=rb.a,
                atom_j=rb.b,
                tag_i=rb.tag_a,
                tag_j=rb.tag_b,
            )
        )

    return BlockGraph(smiles=smiles, blocks=blocks, edges=edges)


# ===================================================================== #
#  Correctness helper: lossless round trip                              #
# ===================================================================== #
def reassemble_canonical(mol: Chem.Mol, bg: BlockGraph) -> str:
    """Rebuild the molecule from blocks + edges and return canonical SMILES.

    Proves the block graph is a lossless partition of the original: combine all
    block atom sets and re-add the cut bonds at the recorded atom pairs. Uses
    the original mol only to copy atom/bond properties.
    """
    rw = Chem.RWMol()
    old2new: dict[int, int] = {}
    for blk in bg.blocks:
        for a in blk.atom_indices:
            old2new[a] = rw.AddAtom(Chem.Atom(mol.GetAtomWithIdx(a)))
    # intra-block bonds (both endpoints in same block)
    block_of = {a: blk.block_id for blk in bg.blocks for a in blk.atom_indices}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in old2new and b in old2new and block_of[a] == block_of[b]:
            rw.AddBond(old2new[a], old2new[b], bond.GetBondType())
    # inter-block bonds (the recorded edges)
    for e in bg.edges:
        orig_bond = mol.GetBondBetweenAtoms(e.atom_i, e.atom_j)
        rw.AddBond(old2new[e.atom_i], old2new[e.atom_j], orig_bond.GetBondType())
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)
