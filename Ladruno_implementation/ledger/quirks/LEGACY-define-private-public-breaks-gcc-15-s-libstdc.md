---
wp: LEGACY
title: "#define private public breaks GCC 15's libstdc++ if it precedes <sstream>"
legacy_seq: 415
---
### `#define private public` breaks GCC 15's libstdc++ if it precedes `<sstream>`

`std::basic_stringbuf::__xfer_bufptrs` is declared `private` and re-declared later;
flipping the keyword makes the second declaration disagree with the first, which
GCC 15 reports as a hard `-Wtemplate-body` ERROR (not a warning). The ADR-97 P2
pre-flight idiom still works — include the standard library and Eigen FIRST, then
`#define private public`, then the project's own headers.
