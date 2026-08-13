# DGGAL IVEA/RTEA oracle generator

This generator treats the BSD-3-Clause DGGAL `dgg` executable as a black box.
It imports no DGGAL implementation code. `profiles.json` records the exact map
from the old registry symbols to concrete DGGAL names, projection, refinement,
root count and identifier profile. The family symbols are intentionally
abstract and generate no cells.

The preferred reference is a clean source build pinned to DGGAL commit
`e16cea7d930e603e09a8310edcd8f58218016e8f`, plus a recorded eC commit. That
DGGAL commit is newer than PyPI 0.0.6 and fixes odd-level RI7H point-to-zone
conversion. Example after building it:

```sh
python3 scripts/oracles/dggal/generate_ivea_rtea_oracles.py \
  --dgg /absolute/path/to/dgg \
  --source-revision e16cea7d930e603e09a8310edcd8f58218016e8f \
  --level 0 --level 1 --lookup-level 0 --lookup-level 1 \
  --output test/oracles/ivea-rtea/dggal-e16cea7
cd test/oracles/ivea-rtea/dggal-e16cea7
shasum -a 256 -c SHA256SUMS
```

The initially committed reconnaissance fixture was made with DGGAL 0.0.6
(`c323c4c444522f16fd6eba0c56ec65714a147c8c`) under Rosetta because its macOS
wheel contains x86_64 runtime libraries. It deliberately seals only level-zero
cells and lookups, which avoid the later odd-level RI7H lookup fix:

```sh
python3 scripts/oracles/dggal/generate_ivea_rtea_oracles.py \
  --runner 'arch -x86_64' \
  --dgg /absolute/path/to/x86_64/venv/bin/dgg \
  --source-revision c323c4c444522f16fd6eba0c56ec65714a147c8c \
  --package-version 0.0.6 \
  --known-limitation 'Predates e16cea7 odd-level RI7H point-to-zone fix; level-zero fixture only' \
  --output test/oracles/ivea-rtea/dggal-0.0.6
```

This reconnaissance manifest's source revision and package version are
caller-asserted and it does not contain an executable hash. Treat the facts as
level-scoped reconnaissance, not release-grade provenance. A replacement made
from the pinned source build must record the eC commit and hash the `dgg`
executable and loaded DGGAL/eC libraries.

`cells.jsonl` contains complete enumerations, packed IDs, WGS84 centroids and
topological vertices, 5x6 chart boundaries, parents, children, and neighbours
with DGGAL's direction labels. Those labels are not asserted to be the package's
canonical CCW order. `lookup.csv` records explicit WGS84 point-to-zone cases,
including orientation and tie probes. Do not use DGGAL as a production runtime
dependency; the sealed files are test facts only.
