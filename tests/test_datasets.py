"""Tests for :mod:`smart_mcts.datasets` (no network required)."""

import re

import pandas as pd
import pytest

from smart_mcts import datasets
from smart_mcts.datasets import _core


def test_data_home_prefers_explicit_arg(tmp_path):
    home = datasets.get_data_home(tmp_path / "here")
    assert home == tmp_path / "here"
    assert home.is_dir()


def test_data_home_uses_env_var(tmp_path, monkeypatch):
    target = tmp_path / "from_env"
    monkeypatch.setenv(_core.DATA_HOME_ENV_VAR, str(target))
    assert datasets.get_data_home() == target
    assert target.is_dir()


def test_data_home_falls_back_to_os_cache(monkeypatch):
    monkeypatch.delenv(_core.DATA_HOME_ENV_VAR, raising=False)
    home = datasets.get_data_home()
    assert "smart_mcts" in str(home)


def test_registry_and_urls_are_consistent():
    assert set(_core.REGISTRY) == set(_core.URLS)
    sha256 = re.compile(r"^sha256:[0-9a-f]{64}$")
    for key, digest in _core.REGISTRY.items():
        assert sha256.match(digest), f"{key}: {digest!r}"
        assert _core.URLS[key].startswith("https://")


def test_pinned_revisions_are_full_sha1():
    for rev in (_core._SMARTCHEMIST_REV, _core._MOSES_REV):
        assert re.match(r"^[0-9a-f]{40}$", rev)


def test_load_moses_rejects_unknown_subset():
    with pytest.raises(ValueError, match="subset must be"):
        datasets.load_moses("validation")


def test_load_moses_subset_filtering(monkeypatch, tmp_path):
    csv = tmp_path / "dataset_v1.csv"
    csv.write_text("SMILES,SPLIT\nCCO,train\nc1ccccc1,test\nCCN,test_scaffolds\n")
    monkeypatch.setattr(datasets.moses, "fetch_file", lambda *a, **k: csv)

    assert datasets.load_moses() == ["CCO", "c1ccccc1", "CCN"]
    assert datasets.load_moses("test") == ["c1ccccc1"]

    frame = datasets.load_moses("train", as_frame=True)
    assert isinstance(frame, pd.DataFrame)
    assert list(frame.columns) == ["SMILES", "SPLIT"]
    assert frame["SMILES"].tolist() == ["CCO"]
