# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [Semantic Versioning](https://semver.org/).

## Unreleased

- Zenodo DOI badge and concept DOI in README and CITATION.cff.

## 0.2.1 (2026-11-26)

Release metadata only. No change to the package's behaviour.

### Added

- Author ORCID in `CITATION.cff` and `.zenodo.json`.
- Guidance files for coding agents: `AGENTS.md` and `CLAUDE.md`.

## 0.2.0 (2026-11-26)

First public release. Earlier versions were developed privately.

### Added

- Public functions `fit_sequences`, `fit_pair` and `fit_wt`.
- Required `pseudocount` and `dispersion` arguments with a single accepted
  value each: `pseudocount="all_zeros"` and `dispersion="fixed"`.
- R reference package in `r/`, which the Python results are checked against.
- Validation against the published tables of Papkou et al. 2023.
- `CITATION.cff`, `.zenodo.json` and `RELEASING.md` for the 0.2.0 release.
- README install line now uses the GitHub tag.

### Fixed

- `docs/validation.md` no longer says the authors' fitness scripts are in the
  Papkou Zenodo record (13837982). They come from the authors' pipeline archive
  on Zenodo, record 8229097 (DOI 10.5281/zenodo.8229097, CC-BY-4.0), and are not
  redistributed here. `lib/fit_model_wt.R` is re-pinned to the Zenodo md5; it
  differs from the earlier local copy only by one trailing space.
