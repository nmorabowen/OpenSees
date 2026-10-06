---
wp: WP-166
title: "WP-166 — Fork artifacts: Linux opensees.so artifact, tag release, installer publish"
pr: "#921"
date: 2026-10-05
status: "in progress"
section: "table"
---
| **WP-166 — Fork artifacts** (apeGmsh program slice nmorabowen/apeGmsh#1488, link F1). `ladruno.yml` job `zone-a-ubuntu` uploads `opensees.so` + `BUILD_INFO.txt` as artifact `opensees-linux-<sha>` (90 days); new tag-triggered `ladruno_release.yml` (`ladruno-v*`) builds Linux and creates a release with `ladruno-opensees-linux-<tag>.tar.gz` (opensees.so, BUILD_INFO.txt, LICENSE from COPYRIGHT); `publish_installer.ps1` attaches the local Inno installer to an existing release, refusing unless `dist\` is a full 5-target build. Recipe in `Ladruno_scripts/README.md`. | CI / release plumbing | **none — no new classes** | `.github/workflows/{ladruno.yml, ladruno_release.yml}`, `Ladruno_scripts/{publish_installer.ps1, make_linux_build_info.sh, README.md}` | **in progress** | [#921](https://github.com/nmorabowen/OpenSees/pull/921) |
