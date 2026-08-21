#!/usr/bin/env python3
"""Generate deterministic rHEALPix and AusPIX black-box oracle vectors."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
from pathlib import Path

import numpy as np
import pyproj

from rhealpixdggs.cell import CELLS0
from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import Ellipsoid, WGS84_A, WGS84_F
from rhealpixdggs.pj_rhealpix import (
    in_rhealpix_image,
    rhealpix_sphere,
    rhealpix_sphere_inverse,
)
from rhealpixdggs.utils import auth_lat


SOURCE_COMMIT = "5929a73b427a33d66a64800051d077ce36bbf901"
SOURCE_URL = "https://github.com/manaakiwhenua/rhealpixdggs-py"
REFERENCE_LICENSE = "MIT (chosen option of LGPL-3.0-or-later OR MIT)"
N_SIDE = 3
N_CHILDREN = N_SIDE**2
BOUNDARY_EDGE_SAMPLES = 5


def make_profiles():
    # max_areal_resolution controls only the reference object's artificial
    # constructor depth cap. The mathematical hierarchy is unbounded.
    sphere = Ellipsoid(R=1.0, radians=True)
    wgs84 = Ellipsoid(a=WGS84_A, f=WGS84_F, lon_0=0, radians=False)
    return {
        "rhealpix_unit_003_00": RHEALPixDGGS(
            sphere, N_side=3, north_square=0, south_square=0,
            max_areal_resolution=1e-24,
        ),
        "auspix_wgs84_003_00": RHEALPixDGGS(
            wgs84, N_side=3, north_square=0, south_square=0,
            max_areal_resolution=1e-6,
        ),
    }


def ordinal0(cell):
    value = CELLS0.index(cell.suid[0])
    for digit in cell.suid[1:]:
        value = value * N_CHILDREN + digit
    return value


def parent_id(cell):
    return None if cell.resolution == 0 else str(cell.rdggs.cell(cell.suid[:-1]))


def children_ids(cell):
    return [str(c) for c in cell.subcells()]


def to_radians(rdggs, point):
    if rdggs.ellipsoid.radians:
        return float(point[0]), float(point[1])
    return tuple(map(float, np.deg2rad(point)))


def unit_xyz(lon_lat_rad):
    lon, lat = lon_lat_rad
    c = math.cos(lat)
    return [c * math.cos(lon), c * math.sin(lon), math.sin(lat)]


def coordinate_views(rdggs, point):
    """Return explicit geodetic and authalic views of an ellipsoid point."""
    geodetic = to_radians(rdggs, point)
    authalic = (
        geodetic[0],
        float(auth_lat(geodetic[1], rdggs.ellipsoid.e, radians=True)),
    )
    return {
        "geodetic_lon_lat_rad": list(geodetic),
        "authalic_lon_lat_rad": list(authalic),
        "unit_direction_authalic": unit_xyz(authalic),
    }


def planar_boundary(cell, n=BOUNDARY_EDGE_SAMPLES):
    """Sample four planar edges CCW from upper-left, without closure."""
    ul = cell.ul_vertex(plane=True)
    delta = cell.width(plane=True) / (n - 1)
    point = np.array(ul, dtype=float)
    result = [tuple(point)]
    # Down, right, up, left is counterclockwise in an x-right/y-up plane.
    for direction in ((0, -1), (1, 0), (0, 1), (-1, 0)):
        origin = point.copy()
        for j in range(1, n):
            result.append(tuple(origin + j * delta * np.array(direction)))
        point = np.array(result[-1])
    result.pop()
    return result


def ellipsoid_boundary(cell, n=BOUNDARY_EDGE_SAMPLES):
    points = [
        cell.rdggs.rhealpix(*p, inverse=True, region=cell.region())
        for p in planar_boundary(cell, n=n)
    ]
    return [coordinate_views(cell.rdggs, p) for p in points]


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path, rows):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def hierarchy_rows(rdggs, max_level=3):
    for level in range(max_level + 1):
        for cell in rdggs.grid(level):
            ul_x, ul_y = cell.ul_vertex(plane=True)
            yield {
                "level": level, "ordinal0": ordinal0(cell), "cell_id": str(cell),
                "parent_id": parent_id(cell) or "",
                "children": " ".join(children_ids(cell)),
                "ul_x": repr(float(ul_x)), "ul_y": repr(float(ul_y)),
                "width": repr(float(cell.width(plane=True))),
            }


def projection_rows():
    transition = math.asin(2 / 3)
    eps = 2.0**-40
    points = [
        (-math.pi, 0.0), (-3 * math.pi / 4, transition),
        (-math.pi / 2, -transition), (-math.pi / 4, transition + eps),
        (0.0, 0.0), (0.0, transition - eps), (0.0, transition),
        (0.0, transition + eps), (math.pi / 4, -transition - eps),
        (math.pi / 2, 0.3), (3 * math.pi / 4, -1.2),
        (math.pi - eps, 0.0), (-2.2, math.pi / 2), (1.1, -math.pi / 2),
    ]
    for north_square in range(4):
        for south_square in range(4):
            for case, (lon, lat) in enumerate(points):
                x, y = rhealpix_sphere(lon, lat, north_square, south_square)
                ilon, ilat = rhealpix_sphere_inverse(x, y, north_square, south_square)
                yield {
                    "north_square": north_square, "south_square": south_square,
                    "case": case, "lon_rad": repr(float(lon)),
                    "lat_rad": repr(float(lat)), "x_authalic_radii": repr(float(x)),
                    "y_authalic_radii": repr(float(y)),
                    "inverse_lon_rad": repr(float(ilon)),
                    "inverse_lat_rad": repr(float(ilat)),
                }


def auspix_projection_rows(rdggs):
    points = [
        (-180.0, 0.0), (-135.0, 70.0), (-90.0, -70.0), (-45.0, 41.0),
        (0.0, 0.0), (45.0, -41.0), (90.0, 20.0), (135.0, -20.0),
        (149.15772, -35.34385), (151.2093, -33.8688),
        (115.8605, -31.9505), (144.9631, -37.8136),
    ]
    for case, (lon, lat) in enumerate(points):
        x, y = rdggs.rhealpix(lon, lat)
        ilon, ilat = rdggs.rhealpix(x, y, inverse=True)
        yield {
            "case": case, "lon_deg": repr(lon), "lat_deg": repr(lat),
            "x_m": repr(float(x)), "y_m": repr(float(y)),
            "inverse_lon_deg": repr(float(ilon)),
            "inverse_lat_deg": repr(float(ilat)),
        }


def validate_projection_against_proj():
    """Cross-check pure-Python spherical projection against MIT-licensed PROJ."""
    max_error = 0.0
    for north_square in range(4):
        for south_square in range(4):
            oracle = pyproj.Proj(
                proj="rhealpix", R=1, north_square=north_square,
                south_square=south_square,
            )
            for lon in np.linspace(-math.pi + 1e-3, math.pi - 1e-3, 23):
                for lat in np.linspace(-math.pi / 2 + 1e-3, math.pi / 2 - 1e-3, 19):
                    x, y = rhealpix_sphere(lon, lat, north_square, south_square)
                    px, py = oracle(lon, lat, radians=True)
                    max_error = max(max_error, math.hypot(x - px, y - py))
    if max_error > 2e-14:
        raise RuntimeError(f"PROJ cross-check failed: max error {max_error}")
    return max_error


def cell_rows(profiles, max_level=2):
    for profile, rdggs in profiles.items():
        for level in range(max_level + 1):
            for cell in rdggs.grid(level):
                nucleus = coordinate_views(rdggs, cell.nucleus(plane=False))
                boundary = ellipsoid_boundary(cell)
                neighbors = cell.neighbors(plane=True)
                yield {
                    "profile": profile, "level": level,
                    "ordinal0": ordinal0(cell), "cell_id": str(cell),
                    "parent_id": parent_id(cell), "children": children_ids(cell),
                    "ellipsoidal_shape": cell.ellipsoidal_shape,
                    "area": float(rdggs.cell_area(level, plane=False)),
                    "area_unit": "steradian" if rdggs.ellipsoid.R_A == 1 else "m2",
                    "nucleus_geodetic_lon_lat_rad": nucleus["geodetic_lon_lat_rad"],
                    "nucleus_authalic_lon_lat_rad": nucleus["authalic_lon_lat_rad"],
                    "nucleus_unit_direction_authalic": nucleus["unit_direction_authalic"],
                    "boundary_winding": "CCW viewed from outside",
                    "boundary_start": "planar upper-left corner",
                    "boundary_geodetic_lon_lat_rad": [
                        p["geodetic_lon_lat_rad"] for p in boundary
                    ],
                    "boundary_authalic_lon_lat_rad": [
                        p["authalic_lon_lat_rad"] for p in boundary
                    ],
                    "boundary_unit_directions_authalic": [
                        p["unit_direction_authalic"] for p in boundary
                    ],
                    "edge_neighbors": {
                        d: str(neighbors[d]) for d in ("up", "right", "down", "left")
                    },
                }


def centroid_rows(profiles):
    for profile, rdggs in profiles.items():
        selected = list(rdggs.grid(0))
        by_shape = {}
        for cell in rdggs.grid(2):
            bucket = by_shape.setdefault(cell.ellipsoidal_shape, [])
            if len(bucket) < 4:
                bucket.append(cell)
        for shape in sorted(by_shape):
            selected.extend(by_shape[shape])
        seen = set()
        for cell in selected:
            if str(cell) in seen:
                continue
            seen.add(str(cell))
            nucleus = coordinate_views(rdggs, cell.nucleus(False))
            centroid = coordinate_views(rdggs, cell.centroid(False))
            yield {
                "profile": profile, "level": cell.resolution,
                "cell_id": str(cell), "ellipsoidal_shape": cell.ellipsoidal_shape,
                "nucleus_geodetic_lon_lat_rad": nucleus["geodetic_lon_lat_rad"],
                "nucleus_authalic_lon_lat_rad": nucleus["authalic_lon_lat_rad"],
                "nucleus_unit_direction_authalic": nucleus["unit_direction_authalic"],
                "coordinate_centroid_geodetic_lon_lat_rad": centroid["geodetic_lon_lat_rad"],
                "coordinate_centroid_authalic_lon_lat_rad": centroid["authalic_lon_lat_rad"],
                "coordinate_centroid_unit_direction_authalic": centroid["unit_direction_authalic"],
            }


def lookup_rows(profiles):
    levels = (0, 1, 2, 3, 5, 10)
    lons = (-math.pi, -2.4, -math.pi / 2, -0.2, 0.0, 0.7, math.pi / 2, 2.8)
    lats = (-1.4, -0.9, -0.4, 0.0, 0.4, 0.9, 1.4)
    for profile, rdggs in profiles.items():
        for level in levels:
            for lon in lons:
                for lat in lats:
                    p = (lon, lat) if rdggs.ellipsoid.radians else tuple(np.rad2deg((lon, lat)))
                    cell = rdggs.cell_from_point(level, p, plane=False)
                    yield {
                        "profile": profile, "level": level,
                        "geodetic_lon_rad": repr(lon), "geodetic_lat_rad": repr(lat),
                        "cell_id": "" if cell is None else str(cell),
                    }


def boundary_tie_rows(profiles):
    """Probe every edge normal and all 2-D quadrants at every corner."""
    for profile, rdggs in profiles.items():
        for cell in [*rdggs.grid(0), *rdggs.grid(1)]:
            x, y = cell.ul_vertex(plane=True)
            w = cell.width(plane=True)
            # Outward normal relative to the source planar square.
            edges = [
                ("top", x + w / 2, y, 0, 1, "up"),
                ("right", x + w, y - w / 2, 1, 0, "right"),
                ("bottom", x + w / 2, y - w, 0, -1, "down"),
                ("left", x, y - w / 2, -1, 0, "left"),
            ]
            corners = [
                ("upper_left", x, y), ("upper_right", x + w, y),
                ("lower_right", x + w, y - w), ("lower_left", x, y - w),
            ]
            # Large enough to cross the reference image's 1e-15 fuzz, tiny
            # relative to even a level-1 cell, and exactly derived from width.
            delta = w * 2.0**-34
            probes = []
            for name, px, py, nx, ny, direction in edges:
                expected = str(cell.neighbor(direction, plane=True))
                for normal_step, position in ((-1, "inside"), (0, "exact"), (1, "outside")):
                    probes.append({
                        "probe_kind": "edge", "location": name,
                        "offset_x_sign": normal_step * nx,
                        "offset_y_sign": normal_step * ny,
                        "normal_position": position,
                        "expected_edge_neighbor": expected,
                        "x": px + normal_step * nx * delta,
                        "y": py + normal_step * ny * delta,
                    })
            for name, px, py in corners:
                # Four genuine 2-D corner quadrants plus the exact corner.
                for sx, sy in ((-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1)):
                    probes.append({
                        "probe_kind": "corner", "location": name,
                        "offset_x_sign": sx, "offset_y_sign": sy,
                        "normal_position": "exact" if (sx, sy) == (0, 0) else "quadrant",
                        "expected_edge_neighbor": "",
                        "x": px + sx * delta, "y": py + sy * delta,
                    })

            R = rdggs.ellipsoid.R_A
            for level in (0, 1, 2, 5):
                for probe in probes:
                    px, py = probe["x"], probe["y"]
                    in_image = in_rhealpix_image(
                        px / R, py / R, rdggs.north_square, rdggs.south_square
                    )
                    row = {
                        "profile": profile, "lookup_level": level,
                        "source_cell": str(cell), "source_level": cell.resolution,
                        "probe_kind": probe["probe_kind"],
                        "location": probe["location"],
                        "offset_x_sign": probe["offset_x_sign"],
                        "offset_y_sign": probe["offset_y_sign"],
                        "normal_position": probe["normal_position"],
                        "expected_edge_neighbor": probe["expected_edge_neighbor"],
                        "plane_x": repr(float(px)), "plane_y": repr(float(py)),
                        "perturbation": repr(float(delta)),
                        "in_projection_image": str(bool(in_image)).lower(),
                        "geodetic_lon_rad": "", "geodetic_lat_rad": "",
                        "authalic_lon_rad": "", "authalic_lat_rad": "",
                        "cell_id": "",
                    }
                    if in_image:
                        lonlat = rdggs.rhealpix(px, py, inverse=True)
                        views = coordinate_views(rdggs, lonlat)
                        geod = views["geodetic_lon_lat_rad"]
                        auth = views["authalic_lon_lat_rad"]
                        found = rdggs.cell_from_point(level, lonlat, plane=False)
                        row.update({
                            "geodetic_lon_rad": repr(geod[0]),
                            "geodetic_lat_rad": repr(geod[1]),
                            "authalic_lon_rad": repr(auth[0]),
                            "authalic_lat_rad": repr(auth[1]),
                            "cell_id": "" if found is None else str(found),
                        })
                    yield row


def vertex_neighbor_rows(rdggs, max_level=2):
    # Black-box geometric incidence. Results are sets, not an order oracle.
    for level in range(max_level + 1):
        cells = list(rdggs.grid(level))
        corners = {
            str(c): np.array([
                p["unit_direction_authalic"] for p in ellipsoid_boundary(c, n=2)
            ])
            for c in cells
        }
        edge = {str(c): set(map(str, c.neighbors(True).values())) for c in cells}
        for cell in cells:
            cid = str(cell)
            neighbors = set(edge[cid])
            for other in cells:
                oid = str(other)
                if oid == cid or oid in neighbors:
                    continue
                distances = np.linalg.norm(
                    corners[cid][:, None, :] - corners[oid][None, :, :], axis=2
                )
                if float(np.min(distances)) <= 2e-10:
                    neighbors.add(oid)
            yield {
                "level": level, "cell_id": cid,
                "vertex_neighbor_ids_sorted": sorted(neighbors),
            }


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    profiles = make_profiles()
    proj_crosscheck_error = validate_projection_against_proj()

    write_csv(output / "hierarchy.csv",
              ("level", "ordinal0", "cell_id", "parent_id", "children", "ul_x", "ul_y", "width"),
              hierarchy_rows(profiles["rhealpix_unit_003_00"]))
    write_csv(output / "projection_unit_all_polar_placements.csv",
              ("north_square", "south_square", "case", "lon_rad", "lat_rad", "x_authalic_radii", "y_authalic_radii", "inverse_lon_rad", "inverse_lat_rad"),
              projection_rows())
    write_csv(output / "projection_auspix_wgs84.csv",
              ("case", "lon_deg", "lat_deg", "x_m", "y_m", "inverse_lon_deg", "inverse_lat_deg"),
              auspix_projection_rows(profiles["auspix_wgs84_003_00"]))
    write_jsonl(output / "cells.jsonl", cell_rows(profiles))
    write_jsonl(output / "centroids.jsonl", centroid_rows(profiles))
    write_csv(output / "lookup.csv",
              ("profile", "level", "geodetic_lon_rad", "geodetic_lat_rad", "cell_id"),
              lookup_rows(profiles))
    write_csv(output / "boundary_ties.csv",
              ("profile", "lookup_level", "source_cell", "source_level",
               "probe_kind", "location", "offset_x_sign", "offset_y_sign",
               "normal_position", "expected_edge_neighbor", "plane_x", "plane_y",
               "perturbation", "in_projection_image", "geodetic_lon_rad",
               "geodetic_lat_rad", "authalic_lon_rad", "authalic_lat_rad", "cell_id"),
              boundary_tie_rows(profiles))
    write_jsonl(output / "vertex_neighbors.jsonl",
                vertex_neighbor_rows(profiles["rhealpix_unit_003_00"]))

    vector_names = (
        "boundary_ties.csv", "cells.jsonl", "centroids.jsonl", "hierarchy.csv",
        "lookup.csv", "projection_auspix_wgs84.csv",
        "projection_unit_all_polar_placements.csv", "vertex_neighbors.jsonl",
    )
    vectors = [output / name for name in vector_names]
    manifest = {
        "schema_version": 2,
        "source": {"url": SOURCE_URL, "commit": SOURCE_COMMIT,
                   "package_version": importlib.metadata.version("rHEALPixDGGS"),
                   "license_elected": REFERENCE_LICENSE},
        "generator_environment": {"python": platform.python_version(),
                                  "numpy": np.__version__,
                                  "scipy": importlib.metadata.version("scipy"),
                                  "pyproj": importlib.metadata.version("pyproj"),
                                  "PROJ": pyproj.proj_version_str},
        "independent_projection_crosscheck": {
            "implementation": "PROJ rHEALPix (MIT)",
            "max_abs_error_authalic_radii": proj_crosscheck_error,
            "tolerance": 2e-14,
        },
        "profiles": {
            "rhealpix_unit_003_00": {"ellipsoid": "unit sphere", "angle_unit": "radian",
                                     "north_square": 0, "south_square": 0, "N_side": 3},
            "auspix_wgs84_003_00": {"ellipsoid": "WGS84", "input_angle_unit": "degree",
                                    "stored_angle_unit": "radian", "lon_0_degrees": 0,
                                    "north_square": 0, "south_square": 0, "N_side": 3},
        },
        "coordinate_schema": {
            "geodetic_lon_lat_rad": "longitude/geodetic latitude on the profile ellipsoid, radians",
            "authalic_lon_lat_rad": "same longitude/authalic latitude on the equal-area sphere, radians",
            "unit_direction_authalic": "Cartesian unit vector formed from authalic longitude/latitude",
            "boundary_winding": "CCW viewed from outside; first point is planar upper-left; implicit closure",
            "boundary_ties": "edge probes cross the outward planar normal; corner probes include all four 2-D quadrants and exact corner",
        },
        "files": {p.name: sha256(p) for p in vectors},
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    checksummed = sorted([*vectors, manifest_path], key=lambda p: p.name)
    (output / "SHA256SUMS").write_text(
        "".join(f"{sha256(p)}  {p.name}\n" for p in checksummed)
    )


if __name__ == "__main__":
    main()
