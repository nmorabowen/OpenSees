# This folder holds a READ-ONLY documentation snapshot of the Lemaitre validation
# tests/helpers. The CANONICAL, runnable copies live in tests/ and Ladruno_scripts/
# (so CI runs them). Tell pytest to ignore this snapshot so it is never collected
# (which would otherwise clash with the canonical tests/ files of the same name).
collect_ignore_glob = ["*.py"]
