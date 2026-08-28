"""Standalone, dataset-scale SmartChemist annotator.

Decoupled from Django and indexed for throughput. The original SmartChemist
iterates all ~41k SMARTS patterns through `HasSubstructMatch` for *every*
molecule (a per-molecule full-library scan via `AnnotatedPattern.objects.all()`).
This module instead:

  1. Loads + compiles the pattern library ONCE (cached to disk).
  2. Screens candidates per molecule with a cheap element/ring descriptor
     prefilter, then an RDKit pattern-fingerprint subset test
     (`AllProbeBitsMatch`), before any expensive substructure match.
  3. Resolves the overshadowing hierarchy with the bond-count bug fixed.

The public `annotate()` output mirrors the original match dicts so it is a
drop-in replacement for `SmartChemist._match_smarts_patterns`.
"""

from __future__ import annotations

import pickle
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors

# Atomic numbers we track for the descriptor prefilter.
_HALOGENS = frozenset((9, 17, 35, 53))


@dataclass
class _Pattern:
    """One compiled, pre-screened annotation pattern."""

    index: int
    name: str
    group: str          # 'functional_group' | 'cyclic' | 'biological'
    smarts: str
    hierarchy: tuple[int, ...]
    query: Chem.Mol
    n_bonds: int
    fp: object                       # pattern fingerprint (for subset screen)
    # descriptor minima the target must meet for this pattern to be possible
    min_heavy: int
    min_rings: int
    min_counts: tuple[int, ...]      # (C, N, O, S, P, halogen, other)


def _atom_counts(mol: Chem.Mol) -> tuple[int, ...]:
    c = n = o = s = p = hal = 0
    for a in mol.GetAtoms():
        z = a.GetAtomicNum()
        if z == 6: c += 1
        elif z == 7: n += 1
        elif z == 8: o += 1
        elif z == 16: s += 1
        elif z == 15: p += 1
        elif z in _HALOGENS: hal += 1
    heavy = mol.GetNumHeavyAtoms()
    other = heavy - c - n - o - s - p - hal
    return (c, n, o, s, p, hal, other)


class SmartChemistAnnotator:
    """In-memory, fingerprint-screened annotator.

    Build once per process (or load from a cached index file), then call
    `annotate(mol)` for each molecule. Safe to pickle the *path* and rebuild
    per worker, or to share the built object via fork.
    """

    def __init__(self, patterns: list[_Pattern]):
        self._patterns = patterns

    # ---- construction ------------------------------------------------- #
    @classmethod
    def from_csvs(
        cls,
        functional_groups: str | Path | None = None,
        cyclic: str | Path | None = None,
        biologicals: str | Path | None = None,
        hierarchy_csv: str | Path | None = None,
    ) -> SmartChemistAnnotator:
        """Compile patterns from the SmartChemist CSV files.

        `cyclic` rows are whole reference scaffolds parsed with MolFromSmiles
        (matched as substructures); the others are SMARTS.
        `hierarchy_csv` (smarts_with_hierarchy.csv) supplies the overshadow
        hierarchy by row index, keyed to `index`.
        """
        hierarchy_by_index: dict[int, tuple[int, ...]] = {}
        if hierarchy_csv:
            hdf = pd.read_csv(hierarchy_csv, skiprows=1)
            for i, row in hdf.iterrows():
                raw = str(row.get("Hierarchy", "")).strip()
                if raw and raw not in ("[]", "nan"):
                    hierarchy_by_index[int(row.iloc[0])] = tuple(
                        int(x) for x in raw.strip("[]").split(",") if x.strip()
                    )

        patterns: list[_Pattern] = []
        idx = 0
        sources = [
            ("functional_group", functional_groups, Chem.MolFromSmarts),
            ("cyclic", cyclic, Chem.MolFromSmiles),
            ("biological", biologicals, Chem.MolFromSmarts),
        ]
        for group, path, parser in sources:
            if not path:
                continue
            df = pd.read_csv(path, skiprows=1)
            for _, row in df.iterrows():
                smarts = str(row["SMARTS"])
                q = parser(smarts)
                if q is None:
                    idx += 1
                    continue
                try:
                    fp = Chem.PatternFingerprint(q)
                except Exception:
                    idx += 1
                    continue
                patterns.append(
                    _Pattern(
                        index=idx,
                        name=str(row["trivialname"]),
                        group=group,
                        smarts=smarts,
                        hierarchy=hierarchy_by_index.get(idx, ()),
                        query=q,
                        n_bonds=q.GetNumBonds(),
                        fp=fp,
                        min_heavy=q.GetNumHeavyAtoms(),
                        min_rings=rdMolDescriptors.CalcNumRings(q)
                        if group != "functional_group" else 0,
                        min_counts=_atom_counts(q),
                    )
                )
                idx += 1
        return cls(patterns)

    def save(self, path: str | Path) -> None:
        # Fingerprints/Mols pickle fine; this caches the expensive compile step.
        with open(path, "wb") as fh:
            pickle.dump(self._patterns, fh, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path: str | Path) -> SmartChemistAnnotator:
        with open(path, "rb") as fh:
            return cls(pickle.load(fh))

    # ---- matching ----------------------------------------------------- #
    def annotate(self, mol: Chem.Mol, remove_overshadowed: bool = False) -> list[dict]:
        heavy = mol.GetNumHeavyAtoms()
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        counts = _atom_counts(mol)
        mol_fp = Chem.PatternFingerprint(mol)

        matches: list[dict] = []
        for p in self._patterns:
            # (1) cheap descriptor prefilter
            if heavy < p.min_heavy or n_rings < p.min_rings:
                continue
            if any(counts[i] < p.min_counts[i] for i in range(len(counts))):
                continue
            # (2) fingerprint subset screen — pattern bits must be in mol bits
            if not DataStructs.AllProbeBitsMatch(p.fp, mol_fp):
                continue
            # (3) the only expensive call, now on << 5% of the library
            for hit in mol.GetSubstructMatches(p.query, useChirality=True):
                matches.append(
                    {
                        "atom_indices": hit,
                        "name": p.name,
                        "group": p.group,
                        "smarts": p.smarts,
                        "n_bonds": p.n_bonds,
                        "hierarchy": p.hierarchy,
                        "index": p.index,
                        "overshadowed": False,
                    }
                )
        self._resolve_overshadowing(matches)
        if remove_overshadowed:
            matches = [m for m in matches if not m["overshadowed"]]
        return matches

    @staticmethod
    def _resolve_overshadowing(matches: list[dict]) -> None:
        """Mark less-specific matches as overshadowed.

        Three rules, mirroring the original intent:
          (a) explicit hierarchy: a match declares which pattern indices it
              supersedes; any such match sharing an atom is overshadowed.
          (b) strict atom-subset: a match whose atoms are a proper subset of
              another's is overshadowed.
          (c) atom-set tie: the match with FEWER bonds is overshadowed
              (BUGFIX: original `elif bonds2 < bonds1` was dead code).
        """
        # (a) hierarchy
        for m in matches:
            if not m["hierarchy"]:
                continue
            superseded = set(m["hierarchy"])
            atoms_m = set(m["atom_indices"])
            for sub in matches:
                if sub["index"] in superseded and atoms_m & set(sub["atom_indices"]):
                    sub["overshadowed"] = True
        # (b) + (c)
        for m in matches:
            if m["overshadowed"]:
                continue
            am = set(m["atom_indices"])
            for sub in matches:
                if sub is m:
                    continue
                asub = set(sub["atom_indices"])
                if am < asub:                       # m strictly inside sub
                    m["overshadowed"] = True
                    break
                if am == asub and m["n_bonds"] < sub["n_bonds"]:
                    m["overshadowed"] = True
                    break


def annotate_supplier(
    annotator: SmartChemistAnnotator, mols: Iterable[Chem.Mol]
) -> list[list[dict]]:
    return [annotator.annotate(m) for m in mols if m is not None]
