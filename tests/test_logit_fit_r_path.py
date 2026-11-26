"""Pin the IRLS path to R's glm.fit where the 1e-8 convergence test is marginal."""

import numpy as np
import pytest

from barseqfit import _fit

# Papkou2023 WT data: reference GCCGATCTC (rows 0-11) and genotype TTGGAAATT
# (rows 12-23); six replicates before then after. Zero counts: none.
SUCCESS = [
    35613, 27972, 30530, 25930, 28702, 22740, 19543, 20476, 22602, 16820, 25243, 14105,
    487, 432, 447, 418, 501, 468, 10962, 10662, 12750, 9169, 13405, 7537,
]  # fmt: skip
FAILURE = [
    69540904, 57257101, 60113245, 53096458, 59036655, 45855134,
    35433789, 34585844, 40123450, 30165781, 44019368, 24432376,
    69576030, 57284641, 60143328, 53121970, 59064856, 45877406,
    35442370, 34595658, 40133302, 30173432, 44031206, 24438944,
]  # fmt: skip
ELAPSED = 3.5166197103891137

# R 4.1.3: glm(cbind(Nreads, Diff) ~ 1 + type * time, binomial) gives
# iterations=3 and se(typetest:time) = 0.0056373166889358848. The final
# deviance ratio is 9.37e-9, only 6% below epsilon, so a deviance formula that
# differs from R's by ~1e-8 absolute takes a fourth iteration.
R_ITERATIONS = 3
R_SE = 0.0056373166889358848


def _inputs():
    time = np.tile(np.repeat([0.0, ELAPSED], 6), 2)
    kind = np.repeat([0.0, 1.0], 12)
    design = np.column_stack([np.ones(24), kind, time, kind * time])
    return (
        np.array(SUCCESS, dtype=float),
        np.array(FAILURE, dtype=float),
        design,
    )


def test_scalar_fit_matches_r_iterations_and_se():
    success, failure, design = _inputs()
    fit = _fit._logit_fit(success, failure, design)
    assert fit.converged
    assert fit.iterations == R_ITERATIONS
    assert fit.se[3] == pytest.approx(R_SE, abs=1e-9)


def test_batch_fit_matches_r_iterations_and_se():
    success, failure, design = _inputs()
    fit = _fit._logit_fit_batch(success[None], failure[None], design[None])
    assert fit.errors == [None]
    assert fit.converged[0]
    assert fit.iterations[0] == R_ITERATIONS
    assert fit.se[0, 3] == pytest.approx(R_SE, abs=1e-9)


# R 4.1.3 quasibinomial refit on the same data:
#   f <- glm(cbind(Nreads, Diff) ~ 1 + type * time, family = quasibinomial(),
#            data = d, control = glm.control(epsilon = 1e-8, maxit = 25))
#   s <- summary(f); s$dispersion; s$coefficients["typetest:time", "Std. Error"]
# gives the values below. R does not clip the dispersion; here phi_raw > 1 so
# the clipped phi applied by barseqfit equals it and the SEs agree.
R_DISPERSION = 9.9555649416275962
R_SE_QUASI = 0.017787109875131293


def test_scalar_dispersion_and_quasi_se_match_r():
    success, failure, design = _inputs()
    fit = _fit._logit_fit(success, failure, design)
    phi = fit.pearson / fit.df_residual
    assert fit.df_residual == 20
    assert phi > 1
    assert phi == pytest.approx(R_DISPERSION, abs=1e-9)
    assert fit.se[3] * np.sqrt(phi) == pytest.approx(R_SE_QUASI, abs=1e-9)


def test_batch_dispersion_and_quasi_se_match_r():
    success, failure, design = _inputs()
    fit = _fit._logit_fit_batch(success[None], failure[None], design[None])
    phi = fit.pearson[0] / fit.df_residual
    assert fit.errors == [None]
    assert fit.df_residual == 20
    assert phi > 1
    assert phi == pytest.approx(R_DISPERSION, abs=1e-9)
    assert fit.se[0, 3] * np.sqrt(phi) == pytest.approx(R_SE_QUASI, abs=1e-9)
