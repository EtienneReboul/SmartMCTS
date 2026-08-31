# Usage

To use Smart MCTS in a project:

```python
import smart_mcts
```

## Datasets

`smart_mcts.datasets` downloads the required corpora with
[`pooch`](https://www.fatiando.org/pooch/), following the same scheme as
[`scikit-fingerprints`](https://github.com/scikit-fingerprints/scikit-fingerprints):
the data lives **outside the repository**, every file is SHA-256 verified against
a pinned registry, and an already-downloaded file is never fetched twice.

### Where the data goes

Resolution order for the data home (first match wins):

1. an explicit `data_dir=` argument to any loader;
2. the `$SMART_MCTS_DATA` environment variable;
3. the per-user OS cache directory — `pooch.os_cache("smart_mcts")`
   (`~/Library/Caches/smart_mcts` on macOS, `~/.cache/smart_mcts` on Linux).

```python
from smart_mcts.datasets import get_data_home

get_data_home()  # default location
get_data_home("/data/smart_mcts")  # or point it anywhere
```

### MOSES

```python
from smart_mcts.datasets import load_moses

smiles = load_moses("train")  # list[str], ~1.58 M molecules
df = load_moses(as_frame=True)  # DataFrame: SMILES, SPLIT
test = load_moses("test_scaffolds", data_dir="/data/smart_mcts", verbose=True)
```

First call downloads `dataset_v1.csv` (~84 MB) once. `subset` is `None`,
`"train"`, `"test"`, or `"test_scaffolds"`.

### SmartChemist SMARTS library

```python
from smart_mcts.datasets import fetch_smartchemist_smarts, load_smartchemist_annotator

paths = fetch_smartchemist_smarts()  # 4 CSV paths + "license" key
annotator = load_smartchemist_annotator()  # compiled + cached SmartChemistAnnotator
```

`load_smartchemist_annotator()` fetches the four CSVs, compiles the ~41 k-pattern
library once (~30 s), and caches the compiled index next to the CSVs
(`<data_home>/smartchemist/annotator_index.pkl`) so later calls are instant. Pass
`rebuild=True` to recompile, `force_update=True` to re-download.

The SMARTS pattern collection is licensed **CC-BY-ND 4.0**. The upstream license
and attribution text is downloaded with the patterns to
`<data_home>/smartchemist/License_for_patterns_here` (also available on its own
via `fetch_smartchemist_license()`). Per that license, `smart-mcts --version`
prints the required attribution notice; see the README's *Acknowledgments*
section for the citation.

Common keyword arguments across all loaders: `data_dir`, `verbose`,
`force_update`.
