# rHEALPix and AusPIX oracle generator

This generator records black-box output from `rHEALPixDGGS` 0.6.0 at commit
`5929a73b427a33d66a64800051d077ce36bbf901`. That revision is explicitly
dual-licensed LGPL-3.0-or-later **or MIT**; this use elects the MIT option. The
Geoscience Australia AusPIX repository is Apache-2.0.

The two profiles are:

- `rhealpix_unit_003_00`: unit sphere, radians, `(north_square,
  south_square)=(0,0)`, `N_side=3`;
- `auspix_wgs84_003_00`: the AusPIX profile, WGS84 geodetic coordinates,
  Greenwich meridian, `(0,0)`, `N_side=3`.

Create and verify the vectors from a clean checkout with:

```sh
python3.11 -m venv .oracle-venv
.oracle-venv/bin/pip install -r scripts/oracles/rhealpix/requirements.txt
.oracle-venv/bin/python scripts/oracles/rhealpix/generate_rhealpix_oracles.py \
  --output test/oracles/rhealpix
shasum -a 256 -c test/oracles/rhealpix/SHA256SUMS
```

Python 3.11 or newer is required by the pinned package. Do not import the
oracle implementation into production code. The generated files are test
fixtures only.

`cells.jsonl` stores 5 samples per planar edge (16 distinct perimeter samples,
implicit closure). This bypasses a reference-package shortcut that returns
only vertices for some shapes: every recorded point is still obtained by the
reference inverse projection. `vertex_neighbors.jsonl` is a set oracle, not
an order oracle; `edge_neighbors` in `cells.jsonl` preserves the reference
direction labels.
