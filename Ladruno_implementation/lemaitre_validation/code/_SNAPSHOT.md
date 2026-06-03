# code/ — documentation snapshot (read-only)

These are **snapshot copies** bundled with the validation report for self-contained reading.
The **canonical, runnable** copies live in the repo tree (so CI runs them):

| snapshot here | canonical location |
|---|---|
| `lemaitre_ref.py` | `tests/_testbed/lemaitre_ref.py` |
| `test_lemaitre_damage.py` | `tests/test_lemaitre_damage.py` |
| `test_lemaitre_validation.py` | `tests/test_lemaitre_validation.py` |
| `test_lemaitre_vs_asdsteel1d.py` | `tests/test_lemaitre_vs_asdsteel1d.py` |
| `test_lemaitre_cyclic.py` | `tests/test_lemaitre_cyclic.py` |
| `test_lemaitre_notched_bar.py` | `tests/test_lemaitre_notched_bar.py` |
| `lemaitre_validation_figures.py` | `Ladruno_scripts/lemaitre_validation_figures.py` |

Run the canonical copies (the snapshot is import-blocked from pytest via `conftest.py`):

```
pytest tests/test_lemaitre_validation.py tests/test_lemaitre_vs_asdsteel1d.py \
       tests/test_lemaitre_cyclic.py -q          # Zone-A
pytest tests/test_lemaitre_notched_bar.py -q     # Zone-B (gmsh)
python Ladruno_scripts/lemaitre_validation_figures.py <dist\bin>   # regenerate figures/
```
