"""Fingerprint the code of every SRC C/C++ file with comments removed.

A comment-only change (a header stamp) leaves every fingerprint unchanged.
Comments are stripped with the quirk lint's clean(), which also blanks string
contents and preprocessor lines, so the preprocessor lines and string literals are
hashed separately from the raw text.

    python fingerprint.py <worktree> <out.json>
"""
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("cqp", root / "ci" / "check_quirk_patterns.py")
cqp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cqp)

STR = re.compile(r'"(?:\\.|[^"\\\n])*"')
out = {}
for p in sorted((root / "SRC").rglob("*")):
    if p.suffix not in (".cpp", ".h", ".hpp", ".cc", ".cxx") or not p.is_file():
        continue
    t = p.read_text(encoding="utf-8", errors="replace")
    code = " ".join("\n".join(cqp.clean(t)).split())            # tokens, whitespace-normalised
    pre = "\n".join(l.strip() for l in t.splitlines() if l.strip().startswith("#"))
    # string literals outside comments: strip // and /* */ comments first
    nocom = re.sub(r"/\*.*?\*/", " ", t, flags=re.S)
    nocom = re.sub(r"//[^\n]*", " ", nocom)
    strs = "\x00".join(STR.findall(nocom))
    h = hashlib.sha256((code + "\x01" + pre + "\x01" + strs).encode("utf-8")).hexdigest()
    out[p.relative_to(root).as_posix()] = h
Path(sys.argv[2]).write_text(json.dumps(out, indent=0), encoding="utf-8")
print(len(out), "files fingerprinted")
