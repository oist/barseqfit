# Releasing 0.2.0

The owner cuts the release. This file does not tag, publish, or change repository visibility.

Version 0.2.0 is the first public version. The public API is `fit_sequences`, `fit_pair` and `fit_wt`, with `pseudocount` and `dispersion` required, each accepting a single value in this release. Default numerical output is unchanged from the frozen baselines.

Version 0.2.1 changes archive metadata only (ORCID, release date, version number). It exists so that the first Zenodo record carries complete metadata; package behaviour is identical to 0.2.0.

Licence for the files in this tree is MIT, copyright 2026 Okinawa Institute of Science and Technology Graduate University. `CITATION.cff` and `.zenodo.json` use that licence, the current version, and author Igors Dubanevics (Okinawa Institute of Science and Technology Graduate University). The author's ORCID (0000-0002-9989-6649) is recorded in both files since 0.2.1. `repository-code` is `https://github.com/oist/barseqfit`, the public home. The README and `CITATION.cff` cite the repository and carry the Zenodo concept DOI (10.5281/zenodo.24548492, all versions). The version DOI is not stored in the tree, because it cannot sit inside the tag that mints it.

Distribution:

- Python package `barseqfit` on PyPI, published by the owner with trusted publishing from `.github/workflows/release.yml`. That workflow runs only when `github.repository` is `oist/barseqfit`.
- The R package under `r/` is ground truth. It is not published to CRAN.
- A GitHub release of the public repository mints the Zenodo DOI when the Zenodo GitHub integration is enabled. `.zenodo.json` supplies the deposit metadata.

## Owner checklist

- [ ] Make the repository public as `oist/barseqfit`, a fresh repository whose history starts with one clean commit.
- [x] Enable the Zenodo GitHub integration on that repository.
- [ ] Configure the PyPI pending trusted publisher for `oist/barseqfit`.
- [ ] Tag `v0.2.0` on the commit that should be the release and create the GitHub release. The tag name without the `v` must equal `0.2.0` in `pyproject.toml`.
- [x] Confirm the Zenodo DOI for that release (0.2.1: 10.5281/zenodo.24548493).
- [x] Add the DOI badge to the README and the DOI to `CITATION.cff` in a following commit.
- [ ] Confirm the PyPI upload from the release workflow (the applicable registry is PyPI; CRAN does not apply).
