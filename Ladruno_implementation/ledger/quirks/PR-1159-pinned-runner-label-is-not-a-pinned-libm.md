---
wp: PR-1159
title: "A pinned runner label is not a pinned libm: `ubuntu-24.04` moved glibc 8.9 -> 8.12 and the bit-exact Concrete3D glibc replay went red"
pr: "#1159"
date: 2026-11-26
---
### A pinned runner label is not a pinned libm: `ubuntu-24.04` moved glibc 8.9 -> 8.12 and the bit-exact Concrete3D glibc replay went red

- **Bites:** the `concrete3d-replay` job ("Concrete3D record replay (glibc …)")
  went red on `ladruno` and every PR from 2026-11-20 on, on the same SHA that
  had passed (issue #1135). Its gate line: "REPLAY ERROR: records written under
  libm=linux-glibc-2.39-0ubuntu8.9 do not replay bit-exactly under
  linux-glibc-2.39-0ubuntu8.12". The kernel blob was `94dbd6a85b…` in the green
  and the red runs alike; only the runner image changed (ubuntu24/20261111.370
  -> 20261116.376).
- **Why:** WP-225 pinned every job to `runs-on: ubuntu-24.04`, which pins the
  distro and the compiler major version, not the image. GitHub rebuilds the
  image under that label every week or so, and a rebuild may take a glibc
  point release. A fixture that is bit-exact on libm (the libm tag
  `linux-glibc-<ldd --version>` is stamped into it, and the guard refuses a
  different one) breaks on the first such rebuild. Installing g++ with apt in a
  plain `ubuntu:24.04` container does not fix it either: `libc6-dev` requires
  the archive's newest `libc6` and upgrades it.
- **Fix/status:** the job runs in `container: buildpack-deps:noble@sha256:…`,
  which ships g++ and git, so there is no apt step and glibc changes only when
  someone edits the digest. The glibc fixture was re-recorded under that image
  (#1159). The re-record procedure for a deliberate digest bump is in the job's
  comment in `ladruno.yml`. Any future libm-bit-exact gate needs the same
  treatment: a digest-pinned container, not a runner label.
