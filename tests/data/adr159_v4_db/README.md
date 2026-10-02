# ADR-159 v4 definitions-stream fixture

A `database File` written by the pre-ADR-159 binary: ladruno `3144e19ba`
(`ladrunoBuild` = `fc75db7f3`, the same `SRC/`), whose contact definitions
stream is version 4 (64 slots; `db.VECs.64.1`). ADR-159 bumped it to 5
(66 slots, `smoothN`/`smoothT` tail-appended) and kept the v4 read.

The model is `_v4_model()` in `tests/test_adr159_mortar_smooth_contact.py`: the
cohesive unit block, `-mortar -augment never -consistanttan -cohesion 2`, no
loads. `test_adr159_restores_a_v4_stream` restores it on the current binary and
requires the analysis to equal a freshly built twin bit for bit.

Regenerate only on a v4 binary, from `tests/`:

    python -c "import test_adr159_mortar_smooth_contact as t; t.write_v4_fixture('data/adr159_v4_db')"

The files are binary (`.gitattributes`); never let git normalise them.
