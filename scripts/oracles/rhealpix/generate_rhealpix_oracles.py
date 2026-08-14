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

from rhealpixdggs.cell import CELLS0
from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import Ellipsoid, WGS84_A, WGS84_F
from rhealpixdggs.pj_rhealpix import rhealpix_sphere, rhealpix_sphere_inverse


SOURCE_COMMIT = "5929a73b427a33d66a64800051d077ce36bbf901"
SOURCE_URL = "https://github.com/manaakiwhenua/rhealpixdggs-py"
REFERENCE_LICENSE = "MIT (chosen option of LGPL-3.0-or-later OR MIT)"
N_SIDE = 3
N_CHILDREN = N_SIDE**2
BOUNDARY_EDGE_SAMPLES = 5


def make_profiles():
    # max_areal_resolution controls only the reference object's artificial
    # constructor depth cap.  The mathematical hierarchy is unbounded.
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


def parse_id(cell_id: str):
    return tuple([cell_id[0], *map(int, cell_id[1:])])


def ordinal0(cell) -> int:
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


def planar_boundary(cell, n=BOUNDARY_EDGE_SAMPLES):
    """Sample all four straight planar edges clockwise, without closure."""
    ul = cell.ul_vertex(plane=True)
    width = cell.width(plane=True)
    delta = width / (n - 1)
    point = np.array(ul, dtype=float)
    result = [tuple(point)]
    for direction in ((1, 0), (0, -1), (-1, 0), (0, 1)):
        origin = point.copy()
        for j in range(1, n):
            result.append(tuple(origin + j * delta * np.array(direction)))
        point = np.array(result[-1])
    result.pop()  # implicit closure
    return result


def ellipsoid_boundary(cell, n=BOUNDARY_EDGE_SAMPLES):
    region = cell.region()
    points = [
        cell.rdggs.rhealpix(*point, inverse=True, region=region)
        for point in planar_boundary(cell, n=n)
    ]
    return [list(to_radians(cell.rdggs, point)) for point in points]


def write_csv(path: Path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")))
            stream.write("\n")


def hierarchy_rows(rdggs, max_level=3):
    for level in range(max_level + 1):
        for cell in rdggs.grid(level):
            ul_x, ul_y = cell.ul_vertex(plane=True)
            yield {
                "level": level,
                "ordinal0": ordinal0(cell),
                "cell_id": str(cell),
                "parent_id": parent_id(cell) or "",
                "children": " ".join(children_ids(cell)),
                "ul_x": repr(float(ul_x)),
                "ul_y": repr(float(ul_y)),
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
        (math.pi - eps, 0.0), (-2.2, math.pi / 2),
        (1.1, -math.pi / 2),
    ]
    for north_square in range(4):
        for south_square in range(4):
            for case, (lon, lat) in enumerate(points):
                x, y = rhealpix_sphere(
                    lon, lat, north_square=north_square,
                    south_square=south_square,
                )
                inverse_lon, inverse_lat = rhealpix_sphere_inverse(
                    x, y, north_square=north_square,
                    south_square=south_square,
                )
                yield {
                    "north_square": north_square,
                    "south_square": south_square,
                    "case": case,
                    "lon_rad": repr(float(lon)),
                    "lat_rad": repr(float(lat)),
                    "x_authalic_radii": repr(float(x)),
                    "y_authalic_radii": repr(float(y)),
                    "inverse_lon_rad": repr(float(inverse_lon)),
                    "inverse_lat_rad": repr(float(inverse_lat)),
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
        inverse_lon, inverse_lat = rdggs.rhealpix(x, y, inverse=True)
        yield {
            "case": case,
            "lon_deg": repr(lon),
            "lat_deg": repr(lat),
            "x_m": repr(float(x)),
            "y_m": repr(float(y)),
            "inverse_lon_deg": repr(float(inverse_lon)),
            "inverse_lat_deg": repr(float(inverse_lat)),
        }


def cell_rows(profiles, max_level=2):
    for profile, rdggs in profiles.items():
        for level in range(max_level + 1):
            for cell in rdggs.grid(level):
                nucleus = to_radians(rdggs, cell.nucleus(plane=False))
                planar_neighbors = cell.neighbors(plane=True)
                yield {
                    "profile": profile,
                    "level": level,
                    "ordinal0": ordinal0(cell),
                    "cell_id": str(cell),
                    "parent_id": parent_id(cell),
                    "children": children_ids(cell),
                    "ellipsoidal_shape": cell.ellipsoidal_shape,
                    "area": float(rdggs.cell_area(level, plane=False)),
                    "nucleus_lon_lat_rad": list(nucleus),
                    "nucleus_unit_direction": unit_xyz(nucleus),
                    "boundary_lon_lat_rad": ellipsoid_boundary(cell),
                    "edge_neighbors": {
                        direction: str(planar_neighbors[direction])
                        for direction in ("up", "right", "down", "left")
                    },
                }


def centroid_rows(profiles):
    # Four deterministic examples of each available ellipsoidal shape at L2,
    # plus all roots.  Skew-quad/dart centroids are numerical integrations in
    # the reference implementation and should be tested with a tolerance.
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
            nucleus = to_radians(rdggs, cell.nucleus(plane=False))
            centroid = to_radians(rdggs, cell.centroid(plane=False))
            yield {
                "profile": profile,
                "level": cell.resolution,
                "cell_id": str(cell),
                "ellipsoidal_shape": cell.ellipsoidal_shape,
                "nucleus_lon_lat_rad": list(nucleus),
                "centroid_lon_lat_rad": list(centroid),
            }


def lookup_rows(profiles):
    levels = (0, 1, 2, 3, 5, 10)
    longitudes = (-math.pi, -2.4, -math.pi / 2, -0.2, 0.0, 0.7, math.pi / 2, 2.8)
    latitudes = (-1.4, -0.9, -0.4, 0.0, 0.4, 0.9, 1.4)
    for profile, rdggs in profiles.items():
        for level in levels:
            for lon in longitudes:
                for lat in latitudes:
                    point = (lon, lat) if rdggs.ellipsoid.radians else tuple(np.rad2deg((lon, lat)))
                    cell = rdggs.cell_from_point(level, point, plane=False)
                    yield {
                        "profile": profile,
                        "level": level,
                        "lon_rad": repr(lon),
                        "lat_rad": repr(lat),
                        "cell_id": "" if cell is None else str(cell),
                    }


def boundary_tie_rows(profiles):
    # Reference ownership on exact boundaries and its two adjacent floating
    # values.  The probes originate in planar space so both root and child
    # seams are covered without reimplementing the classifier.
    for profile, rdggs in profiles.items():
        roots_and_level1 = [*rdggs.grid(0), *rdggs.grid(1)]
        probes = []
        for cell in roots_and_level1:
            ul_x, ul_y = cell.ul_vertex(plane=True)
            width = cell.width(plane=True)
            probes.extend([
                (str(cell), "upper_left", ul_x, ul_y),
                (str(cell), "upper_mid", ul_x + width / 2, ul_y),
                (str(cell), "left_mid", ul_x, ul_y - width / 2),
                (str(cell), "lower_right", ul_x + width, ul_y - width),
            ])
        for level in (0, 1, 2, 5):
            for source_id, location, x, y in probes:
                for side, px in (
                    ("below", np.nextafter(x, -math.inf)),
                    ("exact", x),
                    ("above", np.nextafter(x, math.inf)),
                ):
                    # A single-x ULP probe complements the exact 2-D corner.
                    lonlat = rdggs.rhealpix(px, y, inverse=True)
                    cell = rdggs.cell_from_point(level, lonlat, plane=False)
                    lon_rad, lat_rad = to_radians(rdggs, lonlat)
                    yield {
                        "profile": profile,
                        "level": level,
                        "source_cell": source_id,
                        "location": location,
                        "x_side": side,
                        "lon_rad": repr(lon_rad),
                        "lat_rad": repr(lat_rad),
                        "cell_id": "" if cell is None else str(cell),
                    }


def vertex_neighbor_rows(rdggs, max_level=2):
    # A black-box geometric incidence oracle.  Results are sorted as sets and
    # intentionally do not choose the production API's rotational order.
    for level in range(max_level + 1):
        cells = list(rdggs.grid(level))
        corners = {
            str(cell): np.array([unit_xyz(p) for p in ellipsoid_boundary(cell, n=2)])
            for cell in cells
        }
        edge = {
            str(cell): set(map(str, cell.neighbors(plane=True).values()))
            for cell in cells
        }
        for cell in cells:
            cell_id = str(cell)
            neighbors = set(edge[cell_id])
            for other in cells:
                other_id = str(other)
                if other_id == cell_id or other_id in neighbors:
                    continue
                distances = np.linalg.norm(
                    corners[cell_id][:, None, :] - corners[other_id][None, :, :], axis=2
                )
                if float(np.min(distances)) <= 2e-10:
                    neighbors.add(other_id)
            yield {
                "level": level,
                "cell_id": cell_id,
                "vertex_neighbor_ids_sorted": sorted(neighbors),
            }


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    profiles = make_profiles()

    write_csv(
        output / "hierarchy.csv",
        ("level", "ordinal0", "cell_id", "parent_id", "children", "ul_x", "ul_y", "width"),
        hierarchy_rows(profiles["rhealpix_unit_003_00"]),
    )
    write_csv(
        output / "projection_unit_all_polar_placements.csv",
        ("north_square", "south_square", "case", "lon_rad", "lat_rad", "x_authalic_radii", "y_authalic_radii", "inverse_lon_rad", "inverse_lat_rad"),
        projection_rows(),
    )
    write_csv(
        output / "projection_auspix_wgs84.csv",
        ("case", "lon_deg", "lat_deg", "x_m", "y_m", "inverse_lon_deg", "inverse_lat_deg"),
        auspix_projection_rows(profiles["auspix_wgs84_003_00"]),
    )
    write_jsonl(output / "cells.jsonl", cell_rows(profiles))
    write_jsonl(output / "centroids.jsonl", centroid_rows(profiles))
    write_csv(
        output / "lookup.csv",
        ("profile", "level", "lon_rad", "lat_rad", "cell_id"),
        lookup_rows(profiles),
    )
    write_csv(
        output / "boundary_ties.csv",
        ("profile", "level", "source_cell", "location", "x_side", "lon_rad", "lat_rad", "cell_id"),
        boundary_tie_rows(profiles),
    )
    write_jsonl(
        output / "vertex_neighbors.jsonl",
        vertex_neighbor_rows(profiles["rhealpix_unit_003_00"]),
    )

    vectors = sorted(p for p in output.iterdir() if p.is_file())
    manifest = {
        "schema_version": 1,
        "source": {
            "url": SOURCE_URL,
            "commit": SOURCE_COMMIT,
            "package_version": importlib.metadata.version("rHEALPixDGGS"),
            "license_elected": REFERENCE_LICENSE,
        },
        "generator_environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": importlib.metadata.version("scipy"),
            "pyproj": importlib.metadata.version("pyproj"),
        },
        "profiles": {
            "rhealpix_unit_003_00": {
                "ellipsoid": "unit sphere", "angle_unit": "radian",
                "north_square": 0, "south_square": 0, "N_side": 3,
            },
            "auspix_wgs84_003_00": {
                "ellipsoid": "WGS84", "input_angle_unit": "degree",
                "stored_angle_unit": "radian", "lon_0_degrees": 0,
                "north_square": 0, "south_square": 0, "N_side": 3,
            },
        },
        "files": {p.name: sha256(p) for p in vectors},
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checksummed = sorted([*vectors, manifest_path], key=lambda p: p.name)
    (output / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in checksummed),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
