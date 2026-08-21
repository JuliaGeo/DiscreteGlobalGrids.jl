# IVEA / RTEA implementation scope

This module implements the common 120-fundamental-triangle slice-and-dice
equal-area projection and the registered rhombic `4R` and `9R` refinements for
both profiles.  It uses only the cited papers and BSD-3-Clause DGGAL tables.
Runtime results are computed analytically; oracle files are test inputs only.

The registered `3H` and `7H` variants are intentionally not exposed yet.  The
papers define the Goldberg tessellations, but the current committed DGGAL
oracle is complete only at level zero and does not determine the registered
deep hierarchy.  Correct implementation still needs, from the pinned
post-0.0.6 DGGAL source build:

1. the parity-dependent RI3H/RI7H atlas coordinate transforms and seam ties;
2. the bijection between the twelve pentagonal roots and ordinary atlas zones;
3. canonical ZIRS enumeration (plus the RI7H `_Z7` alternate codec);
4. primary-parent versus all-covering-parent mappings at several deep levels;
5. inverse-projected polygon vertices and point-to-zone probes, especially
   odd-level aperture-7 cases fixed by DGGAL commit `e16cea7`.

Implementing an arbitrary Goldberg grid without those facts would be
mathematically plausible but would not establish compatibility with the old
registry's DGGAL families, so this module does not pretend that it does.
