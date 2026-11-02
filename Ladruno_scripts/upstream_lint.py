"""Upstream PR lint for the jaabell/OpenSees campaign (rules 1 and 6).

Checks the lines a branch adds over its base, and its commit messages, for
fork-internal references, AI traces and tone violations.

    python upstream_lint.py <repo> <base> <head> [--body PR_BODY.md]

Exit 1 when anything is flagged. Ladruno class and command names are allowed
(owner decision 2026-10-30); only internal references are flagged.
"""
import argparse
import re
import subprocess
import sys

# (label, regex, applies_to) ; applies_to: "code", "text" (commits/body), "all"
RULES = [
    ("ADR reference", r"\bADR[- ]?\d+", "all"),
    ("WP reference", r"\bWP[- ]?\d+", "all"),
    ("finding id F##", r"\bF\d{1,3}\b", "text"),
    ("finding id F## in comment", r"(//|/\*|^\s*\*).*(\(F\d{1,3}\)|[Ff]inding\s+F\d{1,3})", "code"),
    ("fork PR/issue number", r"(?<![\w/&])#\d{3,4}\b|nmorabowen/OpenSees", "all"),
    ("apeGmsh / TIMs reference", r"apeGmsh|\bTIMs?\b", "all"),
    ("ledger / quirk reference", r"\bLEDGER\w*|\bquirk\b|LEDGER_", "all"),
    ("Ladruno marker comment", r"//\s*Ladruno\b(?!\w)", "code"),
    ("Guppi", r"Guppi", "all"),
    ("AI trace", r"Claude|Anthropic|Generated with|noreply@anthropic", "all"),
    ("fork-internal path", r"Ladruno_(implementation|internal|scripts)|\.claude/|AGENTS\.md|CLAUDE\.md", "all"),
    ("em dash", r"—", "all"),
    ("first/second person", r"\b(we|We|our|Our|you|You|your|Your)\b", "all"),
    ("bold lead-in", r"^\s*[-*]?\s*\*\*[^*]+\*\*\s*[:.—-]", "text"),
]

SKIP_FILE = re.compile(r"\.(png|jpg|h5|hdf5|bin|gz)$", re.I)


def run(repo, *args):
    return subprocess.run(["git", "-C", repo, *args], capture_output=True,
                          text=True, encoding="utf-8", errors="replace").stdout


def check(text, where, kind, hits):
    for i, line in enumerate(text.splitlines(), 1):
        for label, rx, scope in RULES:
            if scope != "all" and scope != kind:
                continue
            if re.search(rx, line):
                hits.append(f"{where}:{i}: [{label}] {line.strip()[:140]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("repo")
    ap.add_argument("base")
    ap.add_argument("head")
    ap.add_argument("--body")
    a = ap.parse_args()
    hits = []

    diff = run(a.repo, "diff", "-U0", f"{a.base}...{a.head}")
    cur = None
    lineno = 0
    for line in diff.splitlines():
        if line.startswith("+++ "):
            cur = line[6:] if line.startswith("+++ b/") else None
        elif line.startswith("@@"):
            m = re.search(r"\+(\d+)", line)
            lineno = int(m.group(1)) if m else 0
        elif line.startswith("+") and cur and not SKIP_FILE.search(cur):
            body = line[1:]
            for label, rx, scope in RULES:
                if scope == "text":
                    continue
                if re.search(rx, body):
                    hits.append(f"{cur}:{lineno}: [{label}] {body.strip()[:140]}")
            lineno += 1

    log = run(a.repo, "log", "--format=%H%n%an <%ae>%n%B%n---END---", f"{a.base}..{a.head}")
    for chunk in log.split("---END---"):
        chunk = chunk.strip()
        if not chunk:
            continue
        sha, author, msg = (chunk.split("\n", 2) + ["", ""])[:3]
        if author != "Nicolas Mora Bowen <nmorabowen@gmail.com>":
            hits.append(f"commit {sha[:9]}: [author] {author}")
        for need in ("Co-authored-by: Patricio Palacios <pxpalacios@miuandes.cl>",
                     "Co-authored-by: Jose A. Abell <jaabell@miuandes.cl>"):
            if need not in msg:
                hits.append(f"commit {sha[:9]}: [missing trailer] {need}")
        check(msg, f"commit {sha[:9]}", "text", hits)

    if a.body:
        with open(a.body, encoding="utf-8") as f:
            check(f.read(), a.body, "text", hits)

    for h in hits:
        print(h)
    print(f"{len(hits)} finding(s)")
    sys.exit(1 if hits else 0)


if __name__ == "__main__":
    main()
