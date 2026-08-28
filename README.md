# Smart MCTS

[![PyPI version](https://img.shields.io/pypi/v/smart-mcts.svg)](https://pypi.org/project/smart-mcts/)
[![PyPI downloads](https://static.pepy.tech/badge/smart-mcts/month)](https://pepy.tech/projects/smart-mcts)

Graph-based Monte Carlo Tree Search that assembles drug-like molecules from chemically meaningful blocks, guided by a PUCT prior learned from real datasets.

* [GitHub](https://github.com/EtienneReboul/smart-mcts/) | [PyPI](https://pypi.org/project/smart-mcts/) | [Documentation](https://EtienneReboul.github.io/smart-mcts/)
* Created by [Etienne Reboul](https://github.com/EtienneReboul) | GitHub [@EtienneReboul](https://github.com/EtienneReboul) | PyPI [@EtienneReboul](https://pypi.org/user/EtienneReboul/)
* MIT License

## Features

* `smart_mcts.smartchemist_annotator` — Django-free, fingerprint-screened
  re-implementation of SmartChemist hierarchical SMARTS annotation (~5x faster
  than a naive full-library scan; results verified identical).
* `smart_mcts.hybrid_joiner` — turns a molecule into a block graph whose nodes
  are annotated motifs and whose edges are R-BRICS-breakable bonds between
  motifs; reassembly is verified lossless.
* Emits tag-annotated transition records for learning the
  `P(next_block | current_block, bond_tag)` prior that guides the PUCT search.

See [handoff.md](handoff.md) for full project context, design decisions, and the
prioritised next steps.

## Installation

```bash
uv add smart-mcts
```

## Usage

```python
import smart_mcts
```

## Documentation

Full documentation is available on
[GitHub Pages](https://EtienneReboul.github.io/smart-mcts/).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, testing, and
documentation instructions.

## Author

Smart MCTS was created in 2026 by Etienne Reboul.

Built with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
