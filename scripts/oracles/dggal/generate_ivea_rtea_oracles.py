#!/usr/bin/env python3
"""Seal IVEA/RTEA facts from the BSD-3-Clause DGGAL CLI.

This script treats ``dgg`` as an external process.  It intentionally imports no
DGGAL implementation code and records only observable identifiers, geometry and
relationships.  The caller must state the exact source revision being run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import shlex
import subprocess
from pathlib import Path


SOURCE_URL = "https://github.com/ecere/dggal"
SOURCE_LICENSE = "BSD-3-Clause"
ORIENTATION = {
    "model": "DGGAL/OGC IVEA 5x6 atlas",
    "first_vertex_authalic_lat_degrees": 58.282525588538994,
    "first_vertex_geodetic_wgs84_lat_degrees": 58.397145907431,
    "first_vertex_lon_degrees_east": 11.2,
    "adjacent_designated_vertex": "due north",
    "cli_geographic_crs": "EPSG:4326 (latitude,longitude argument order)",
    "planar_crs": "DGGAL 5x6",
}

# This is also the explicit mapping from every old registry name.  Family names
# are abstract and therefore deliberately have no dggal_cli value.
PROFILES = {
    "IVEA_family": {"concrete": False, "projection": "IVEA"},
    "IVEA4R": {"concrete": True, "dggal_cli": "ivea4r", "projection": "IVEA", "refinement": "RI4R", "aperture": 4, "cell": "rhomb", "root_count": 10, "index": "DGGAL RI4R ZIRS"},
    "IVEA9R": {"concrete": True, "dggal_cli": "ivea9r", "projection": "IVEA", "refinement": "RI9R", "aperture": 9, "cell": "rhomb", "root_count": 10, "index": "OGC/DGGAL IVEA9R ZIRS"},
    "IVEA3H": {"concrete": True, "dggal_cli": "ivea3h", "projection": "IVEA", "refinement": "RI3H", "aperture": 3, "cell": "hexagon/pentagon", "root_count": 12, "index": "OGC/DGGAL IVEA3H ZIRS"},
    "IVEA7H": {"concrete": True, "dggal_cli": "ivea7h", "projection": "IVEA", "refinement": "RI7H", "aperture": 7, "cell": "hexagon/pentagon", "root_count": 12, "index": "OGC/DGGAL IVEA7H ZIRS"},
    "IVEA7H_Z7": {"concrete": True, "dggal_cli": "ivea7h_z7", "projection": "IVEA", "refinement": "RI7H", "aperture": 7, "cell": "hexagon/pentagon", "root_count": 12, "index": "Z7 alternate encoding"},
    "RTEA_family": {"concrete": False, "projection": "DGGAL RT(S)EA"},
    "RTEA4R": {"concrete": True, "dggal_cli": "rtea4r", "projection": "DGGAL RT(S)EA", "refinement": "RI4R", "aperture": 4, "cell": "rhomb", "root_count": 10, "index": "DGGAL RI4R ZIRS"},
    "RTEA9R": {"concrete": True, "dggal_cli": "rtea9r", "projection": "DGGAL RT(S)EA", "refinement": "RI9R", "aperture": 9, "cell": "rhomb", "root_count": 10, "index": "DGGAL RI9R ZIRS"},
    "RTEA3H": {"concrete": True, "dggal_cli": "rtea3h", "projection": "DGGAL RT(S)EA", "refinement": "RI3H", "aperture": 3, "cell": "hexagon/pentagon", "root_count": 12, "index": "DGGAL RI3H ZIRS"},
    "RTEA7H": {"concrete": True, "dggal_cli": "rtea7h", "projection": "DGGAL RT(S)EA", "refinement": "RI7H", "aperture": 7, "cell": "hexagon/pentagon", "root_count": 12, "index": "DGGAL RI7H ZIRS"},
    "RTEA7H_Z7": {"concrete": True, "dggal_cli": "rtea7h_z7", "projection": "DGGAL RT(S)EA", "refinement": "RI7H", "aperture": 7, "cell": "hexagon/pentagon", "root_count": 12, "index": "Z7 alternate encoding"},
}

LOOKUP_POINTS = [
    ("north_pole", 90.0, 0.0),
    ("south_pole", -90.0, 0.0),
    ("origin", 0.0, 0.0),
    ("equator_east", 0.0, 90.0),
    ("equator_antimeridian", 0.0, 180.0),
    ("first_vertex_wgs84", 58.397145907431, 11.2),
    ("first_vertex_minus", 58.397145907431 - 1e-10, 11.2 - 1e-10),
    ("first_vertex_plus", 58.397145907431 + 1e-10, 11.2 + 1e-10),
    ("sydney", -33.8688, 151.2093),
    ("new_york", 40.7128, -74.0060),
]


def command(args, *tail):
    return [*shlex.split(args.runner), str(args.dgg), *map(str, tail)]


def run(args, *tail):
    result = subprocess.run(
        command(args, *tail), check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    return result.stdout


def one(pattern, text, description, flags=0):
    match = re.search(pattern, text, flags)
    if not match:
        raise ValueError(f"could not parse {description}")
    return match


def section(text, heading):
    match = re.search(rf"^{re.escape(heading)}(?: \(\d+\))?:\n((?:   .*\n)+)", text, re.M)
    if not match:
        return []
    return [line.strip() for line in match.group(1).splitlines()]


def relation(line):
    match = one(r"^(\S+?)(?: \((.*)\))?$", line, "relationship")
    return {"id": match.group(1), "annotation": match.group(2)}


def parse_info(text):
    level = one(r"^Level (\d+) zone \((\d+) edges(?:, ([^)]+))?\)$", text, "level", re.M)
    centre = one(r"^WGS84 Centroid \(lat, lon\): ([^,]+), (.+)$", text, "centroid", re.M)
    area = one(r"^([0-9.eE+-]+) m² \(", text, "area", re.M)
    integer = one(r"^64-bit integer ID: (\d+) \((0x[0-9A-Fa-f]+)\)$", text, "integer ID", re.M)
    text_id = one(r"^Textual Zone ID: (\S+)$", text, "text ID", re.M).group(1)
    neighbours = []
    for line in section(text, "Neighbors"):
        match = one(r"^\(direction (\d+)\): (\S+)$", line, "neighbour")
        neighbours.append({"direction": int(match.group(1)), "id": match.group(2)})
    vertices = []
    vertex_match = re.search(r"^\[EPSG:4326\] Vertices \(\d+\):\n((?:   .*\n?)+)", text, re.M)
    if vertex_match:
        for line in vertex_match.group(1).splitlines():
            lat, lon = map(float, line.strip().split(", "))
            vertices.append([lon, lat])
    parents = [] if "\nNo parent\n" in text else [relation(x) for x in section(text, "Parent")]
    if not parents and "\nNo parent\n" not in text:
        parents = [relation(x) for x in section(text, "Parents")]
    return {
        "id": text_id,
        "integer_id": int(integer.group(1)),
        "integer_id_hex": integer.group(2).lower(),
        "level": int(level.group(1)),
        "edge_count": int(level.group(2)),
        "zone_annotation": level.group(3),
        "area_wgs84_m2": float(area.group(1)),
        "centroid_wgs84_lon_lat_deg": [float(centre.group(2)), float(centre.group(1))],
        "vertices_wgs84_lon_lat_deg": vertices,
        "parents": parents,
        "children": [relation(x) for x in section(text, "Children")],
        "neighbors": neighbours,
    }


def chart_boundary(args, cli_name, zone_id):
    feature = json.loads(run(args, cli_name, "-crs", "5x6", "geom", zone_id))
    return feature["geometry"]["coordinates"][0]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_jsonl(path, values):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for value in values:
            stream.write(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def concrete_profiles(selected):
    names = selected or [name for name, p in PROFILES.items() if p["concrete"]]
    unknown = sorted(set(names) - set(PROFILES))
    abstract = sorted(name for name in names if name in PROFILES and not PROFILES[name]["concrete"])
    if unknown or abstract:
        raise ValueError(f"unknown profiles={unknown}; abstract profiles={abstract}")
    return names


def generate(args):
    args.output.mkdir(parents=True, exist_ok=True)
    profile_names = concrete_profiles(args.profile)
    write_json(args.output / "profiles.json", {name: PROFILES[name] for name in PROFILES})

    cells = []
    for registry_name in profile_names:
        profile = PROFILES[registry_name]
        cli_name = profile["dggal_cli"]
        for requested_level in args.level:
            ids = json.loads(run(args, cli_name, "list", requested_level))
            expected = 10 * profile["aperture"] ** requested_level
            if profile["cell"] == "hexagon/pentagon":
                expected += 2
            if len(ids) != expected:
                raise RuntimeError(f"{registry_name} level {requested_level}: {len(ids)} != {expected}")
            for ordinal, zone_id in enumerate(ids):
                record = parse_info(run(args, cli_name, "info", zone_id))
                if record["id"] != zone_id or record["level"] != requested_level:
                    raise RuntimeError(f"DGGAL identity mismatch for {registry_name}/{zone_id}")
                record.update({
                    "registry_name": registry_name,
                    "dggal_cli_name": cli_name,
                    "enumeration_ordinal": ordinal,
                    "boundary_5x6_xy": chart_boundary(args, cli_name, zone_id),
                })
                cells.append(record)

    # Neighbour IDs must exist in each complete enumerated level.  Direction
    # order is deliberately retained as a labelled oracle, not called CCW.
    by_grid_level = {}
    for record in cells:
        by_grid_level.setdefault((record["registry_name"], record["level"]), set()).add(record["id"])
    for record in cells:
        ids = by_grid_level[(record["registry_name"], record["level"])]
        missing = [n["id"] for n in record["neighbors"] if n["id"] not in ids]
        if missing:
            raise RuntimeError(f"missing same-level neighbours for {record['registry_name']}/{record['id']}: {missing}")
    write_jsonl(args.output / "cells.jsonl", cells)

    lookup_path = args.output / "lookup.csv"
    with lookup_path.open("w", encoding="utf-8", newline="") as stream:
        fields = ["registry_name", "dggal_cli_name", "level", "case", "lat_wgs84_deg", "lon_wgs84_deg", "zone_id"]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for registry_name in profile_names:
            cli_name = PROFILES[registry_name]["dggal_cli"]
            for level in args.lookup_level:
                for case, lat, lon in LOOKUP_POINTS:
                    text = run(args, cli_name, "zone", f"{lat:.15g},{lon:.15g}", level)
                    zone_id = one(r"^Textual Zone ID: (\S+)$", text, "lookup ID", re.M).group(1)
                    writer.writerow({
                        "registry_name": registry_name, "dggal_cli_name": cli_name,
                        "level": level, "case": case,
                        "lat_wgs84_deg": repr(lat), "lon_wgs84_deg": repr(lon),
                        "zone_id": zone_id,
                    })

    files = ["cells.jsonl", "lookup.csv", "profiles.json"]
    manifest = {
        "schema_version": 1,
        "purpose": (
            "level-scoped black-box IVEA/RTEA reconnaissance facts; only levels "
            "listed in complete_cell_levels and lookup_levels are implementation-gating"
        ),
        "source": {
            "name": "DGGAL", "url": SOURCE_URL, "license": SOURCE_LICENSE,
            "revision": args.source_revision, "package_version": args.package_version,
            "known_limitations": args.known_limitation,
        },
        "generator_environment": {
            "python": platform.python_version(), "platform": platform.platform(),
            "command_prefix": args.runner,
        },
        "orientation": ORIENTATION,
        "profiles": profile_names,
        "complete_cell_levels": args.level,
        "lookup_levels": args.lookup_level,
        "neighbor_order_semantics": "DGGAL direction labels; not asserted to be CCW",
        "files": {name: sha256(args.output / name) for name in files},
    }
    write_json(args.output / "manifest.json", manifest)
    files.append("manifest.json")
    with (args.output / "SHA256SUMS").open("w", encoding="ascii", newline="\n") as stream:
        for name in files:
            stream.write(f"{sha256(args.output / name)}  {name}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dgg", type=Path, required=True, help="path to pinned dgg executable")
    parser.add_argument("--runner", default="", help="optional shell-split prefix, e.g. 'arch -x86_64'")
    parser.add_argument("--source-revision", required=True, help="exact DGGAL git commit")
    parser.add_argument("--package-version", default=None)
    parser.add_argument("--known-limitation", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", action="append", choices=sorted(PROFILES))
    parser.add_argument("--level", action="append", type=int, default=[])
    parser.add_argument("--lookup-level", action="append", type=int, default=[])
    args = parser.parse_args()
    args.level = sorted(set(args.level or [0]))
    args.lookup_level = sorted(set(args.lookup_level or [0]))
    generate(args)


if __name__ == "__main__":
    main()
