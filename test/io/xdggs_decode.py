"""Open a DGGS store the way an xdggs user does and report what came back.

Run by ``test/io/xdggs_python.jl`` under the interpreter ``DGG_XDGGS_PYTHON``
names. Each fact is printed as ``key=value`` on its own line; a failure at any
step prints ``error=<type>: <message>`` and exits non-zero.
"""

import sys

import xarray as xr
import xdggs

path = sys.argv[1]
try:
    ds = xr.open_dataset(path, engine="zarr")
    decoded = xdggs.decode(ds)
    info = decoded.dggs.grid_info
    centers = decoded.dggs.cell_centers()
    print(f"grid_name={info.to_dict()['grid_name']}")
    print(f"level={info.level}")
    print(f"indexing_scheme={getattr(info, 'indexing_scheme', '')}")
    print(f"ncells={decoded.sizes['cell_ids']}")
    print(f"dtype={decoded['cell_ids'].dtype}")
    print(f"first_id={int(decoded['cell_ids'][0])}")
    print(f"first_lon={float(centers['longitude'][0])!r}")
    print(f"first_lat={float(centers['latitude'][0])!r}")
    print(f"data_vars={','.join(sorted(map(str, decoded.data_vars)))}")
except Exception as exc:  # noqa: BLE001 - the Julia side reports it verbatim
    print(f"error={type(exc).__name__}: {exc}")
    sys.exit(1)
