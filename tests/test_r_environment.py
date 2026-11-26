import json
from pathlib import Path


def test_r_reference_environment_is_pinned_to_paper_r_version():
    lockfile = Path(__file__).parents[1] / "r" / "renv.lock"
    lock = json.loads(lockfile.read_text())

    assert lock["R"]["Version"] == "4.1.3"
    assert lock["Packages"]["renv"]["Version"] == "1.0.11"
    assert {"renv", "testthat"} <= set(lock["Packages"])
