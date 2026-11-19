---
wp: WP-239
title: "An element that caches another Element* dangles after remove element"
date: 2026-11-19
---
### `remove element` + re-create with the same tag: a cached `Element*` in another element dangles

`Domain::removeElement` deletes nothing itself, but the `remove element`
command deletes the object, and nothing tells the elements that hold a pointer
to it. `LadrunoMembrane` resolved its `-frames` once (in `update()`, called
from `Domain::addElement`) and cached `Element*`; a hinge step that removes each
beam and re-adds it on twin nodes (same tag) left it writing slab loads through
freed memory. It gave the right totals only because glibc handed the freed
block straight back to the new beam of the same size — luck, not design (YOO
RevB NLRHA, 13 292 beams, 2026-11-19). Rules: hold another element by **tag**
and look it up through `getDomain()->getElement(tag)` where it is used (a map
lookup — cheap next to assembly); never keep the pointer across calls. And when
the referenced element may move onto coincident twin nodes, match its ends to
your own nodes by coordinate as well as by tag, or a database `restore` (which
re-resolves against the current topology) builds a different model from the
live one. Fixed for `LadrunoMembrane` in WP-239.
