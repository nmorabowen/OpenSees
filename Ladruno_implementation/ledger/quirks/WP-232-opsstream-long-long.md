---
wp: WP-232
title: "OPS_Stream << long long is ambiguous"
date: 2026-11-13
---
### `OPS_Stream << long long` is ambiguous

- **Bites:** `opserr << (long long)id` does not compile (GCC: ambiguous overload between
  the int and double inserters). 64-bit particle ids hit it.
- **Workaround:** `opserr << std::to_string(id).c_str()`.
- Learned 2026-11-13 (WP-232c).
