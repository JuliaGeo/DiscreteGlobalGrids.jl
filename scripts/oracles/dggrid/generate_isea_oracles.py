#!/usr/bin/env python3
"""Generate sealed ISEA oracle vectors by treating DGGRID as a black box.

No DGGRID source is imported, copied, or consulted by this script.  The caller
supplies a DGGRID executable; this program renders a public CLI meta file,
executes it, validates the resulting record counts, and records the raw output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


PINNED_DGGRID_COMMIT = "04bf5ed1372b174b9349faca8c265d112f6d8587"
PINNED_DGGRID_VERSION = "9.0b"
STANDARD_VERTEX_LATITUDE = 58.282525588538995


@dataclass(frozen=True)
class GridSpec:
    name: str
    max_level: int
    address_type: str
    hierarchy: str
    neighbors: bool
    index_parents: bool

    def count(self, level: int) -> int:
        if self.name == "ISEA3H":
            return 10 * 3**level + 2
        if self.name == "ISEA4H":
            return 10 * 4**level + 2
        if self.name == "ISEA4T":
            return 20 * 4**level
        raise AssertionError(self.name)


GRID_SPECS = (
    GridSpec("ISEA3H", 5, "HIERNDX", "Z3", True, True),
    GridSpec("ISEA4H", 4, "HIERNDX", "ZORDER", True, True),
    GridSpec("ISEA4T", 4, "SEQNUM", "ZORDER", False, False),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def render(template: str, values: dict[str, str]) -> str:
    result = template
    for key, value in values.items():
        result = result.replace(f"@{key}@", value)
    leftovers = sorted(set(re.findall(r"@[A-Z_]+@", result)))
    if leftovers:
        raise RuntimeError(f"unexpanded template fields: {leftovers}")
    return result


def invoke(
    executable: Path,
    template: str,
    spec: GridSpec,
    level: int,
    run_dir: Path,
    *,
    hierarchy_form: str,
    geometry: bool,
    vertex_longitude: float = 11.25,
    vertex_latitude: float = STANDARD_VERTEX_LATITUDE,
    vertex_azimuth: float = 0.0,
) -> tuple[str, str]:
    prefix = run_dir / "oracle"
    values = {
        "DGGS_TYPE": spec.name,
        "RESOLUTION": str(level),
        "VERTEX_LONGITUDE": format(vertex_longitude, ".15f"),
        "VERTEX_LATITUDE": format(vertex_latitude, ".15f"),
        "VERTEX_AZIMUTH": format(vertex_azimuth, ".15f"),
        "ADDRESS_TYPE": spec.address_type,
        "HIER_SYSTEM": spec.hierarchy,
        "HIER_FORM": hierarchy_form,
        "CELL_OUTPUT_TYPE": "GEOJSON" if geometry else "NONE",
        "CELL_FILE": str(prefix) + "-boundaries",
        "POINT_FILE": str(prefix) + "-centers",
        "NEIGHBOR_OUTPUT_TYPE": "TEXT" if geometry and spec.neighbors else "NONE",
        "NEIGHBOR_FILE": str(prefix) + "-neighbors",
        "INDEX_OUTPUT_TYPE": (
            "TEXT" if geometry and spec.index_parents and level > 0 else "NONE"
        ),
        "PARENT_FILE": str(prefix) + "-parents",
        "CHILDREN_FILE": str(prefix) + "-children",
    }
    meta = run_dir / "run.meta"
    meta.write_text(render(template, values), encoding="utf-8")
    process = subprocess.run(
        [str(executable), str(meta)],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        raise RuntimeError(
            f"DGGRID failed for {spec.name} level {level}:\n"
            f"{process.stdout}\n{process.stderr}"
        )
    version_marker = f"DGGRID version {PINNED_DGGRID_VERSION}"
    if version_marker not in process.stdout:
        raise RuntimeError(
            f"expected {version_marker!r} in DGGRID output; refusing an unpinned version"
        )
    expected_marker = f"generated {spec.count(level)} cells"
    if expected_marker not in process.stdout.replace(",", ""):
        raise RuntimeError(
            f"missing {expected_marker!r} for {spec.name} level {level}"
        )
    return process.stdout, process.stderr


def copy_checked(source: Path, destination: Path, expected_lines: int | None) -> dict:
    if not source.is_file():
        raise RuntimeError(f"DGGRID did not create {source}")
    if expected_lines is not None:
        with source.open(encoding="utf-8") as stream:
            lines = sum(1 for _ in stream)
        if lines != expected_lines:
            raise RuntimeError(
                f"{source.name}: expected {expected_lines} lines, found {lines}"
            )
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return {
        "path": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
    }


def generate(executable: Path, output_root: Path) -> None:
    template_path = Path(__file__).with_name("isea.meta.in")
    template = template_path.read_text(encoding="utf-8")

    for spec in GRID_SPECS:
        destination_dir = output_root / spec.name / "dggrid-9.0b"
        destination_dir.mkdir(parents=True, exist_ok=True)
        files: list[dict] = []
        runs: list[dict] = []
        for level in range(spec.max_level + 1):
            expected = spec.count(level)
            with tempfile.TemporaryDirectory(prefix="dgg-oracle-") as temporary:
                run_dir = Path(temporary)
                stdout, stderr = invoke(
                    executable,
                    template,
                    spec,
                    level,
                    run_dir,
                    hierarchy_form="DIGIT_STRING",
                    geometry=True,
                )
                prefix = run_dir / "oracle"
                files.append(
                    copy_checked(
                        Path(str(prefix) + "-centers.txt"),
                        destination_dir / f"level-{level:02d}-centers.txt",
                        expected,
                    )
                )
                boundary_source = Path(str(prefix) + "-boundaries.geojson")
                with boundary_source.open(encoding="utf-8") as stream:
                    feature_count = len(json.load(stream)["features"])
                if feature_count != expected:
                    raise RuntimeError(
                        f"{spec.name} level {level}: expected {expected} boundary features, "
                        f"found {feature_count}"
                    )
                files.append(
                    copy_checked(
                        boundary_source,
                        destination_dir / f"level-{level:02d}-boundaries.geojson",
                        None,
                    )
                )
                if spec.neighbors:
                    files.append(
                        copy_checked(
                            Path(str(prefix) + "-neighbors.nbr"),
                            destination_dir / f"level-{level:02d}-neighbors.txt",
                            expected,
                        )
                    )
                if spec.index_parents and level > 0:
                    files.append(
                        copy_checked(
                            Path(str(prefix) + "-parents.ndxPrt"),
                            destination_dir / f"level-{level:02d}-parents.txt",
                            expected,
                        )
                    )
                runs.append(
                    {
                        "level": level,
                        "cell_count": expected,
                        "stderr": stderr.strip(),
                        "stdout_completion": expected_marker(stdout),
                    }
                )

            if spec.address_type == "HIERNDX":
                with tempfile.TemporaryDirectory(prefix="dgg-oracle-codec-") as temporary:
                    run_dir = Path(temporary)
                    invoke(
                        executable,
                        template,
                        spec,
                        level,
                        run_dir,
                        hierarchy_form="INT64",
                        geometry=False,
                    )
                    files.append(
                        copy_checked(
                            run_dir / "oracle-centers.txt",
                            destination_dir / f"level-{level:02d}-{spec.hierarchy.lower()}-int64-centers.txt",
                            expected,
                        )
                    )

        manifest = {
            "schema": 1,
            "source": "DGGRID CLI black-box output",
            "dggrid_version": PINNED_DGGRID_VERSION,
            "dggrid_commit_used_to_build_local_executable": PINNED_DGGRID_COMMIT,
            "dggrid_executable_sha256": sha256(executable),
            "implementation_reference": False,
            "grid": spec.name,
            "orientation": {
                "vertex_0_longitude_degrees": 11.25,
                "vertex_0_latitude_degrees": STANDARD_VERTEX_LATITUDE,
                "vertex_0_azimuth_degrees": 0.0,
            },
            "sphere_radius": 1.0,
            "precision_digits_after_decimal": 15,
            "address_type": spec.address_type,
            "hierarchy": spec.hierarchy if spec.address_type == "HIERNDX" else None,
            "levels": runs,
            "files": sorted(files, key=lambda item: item["path"]),
        }
        manifest_path = destination_dir / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    generate_ogc_orientation(executable, output_root, template)


def generate_ogc_orientation(executable: Path, output_root: Path, template: str) -> None:
    """Generate a small spherical corpus for the OGC Annex-B orientation.

    Annex B gives 58.397145907431 degrees as a *WGS84 geodetic* latitude.
    Its authalic latitude is atan(phi), the spherical latitude below.  These
    files isolate the 0.05-degree longitude rotation from the package's
    DGGRID-standard 11.25-degree frame; they do not claim to test the OGC ZIRS.
    """

    spec = GRID_SPECS[0]
    longitude = 11.20
    latitude = STANDARD_VERTEX_LATITUDE
    destination_dir = output_root / spec.name / "ogc-annex-b-orientation-dggrid-9.0b"
    destination_dir.mkdir(parents=True, exist_ok=True)
    files: list[dict] = []
    runs: list[dict] = []
    for level in range(4):
        expected = spec.count(level)
        with tempfile.TemporaryDirectory(prefix="dgg-oracle-ogc-orientation-") as temporary:
            run_dir = Path(temporary)
            stdout, stderr = invoke(
                executable,
                template,
                spec,
                level,
                run_dir,
                hierarchy_form="DIGIT_STRING",
                geometry=True,
                vertex_longitude=longitude,
                vertex_latitude=latitude,
            )
            prefix = run_dir / "oracle"
            files.append(
                copy_checked(
                    Path(str(prefix) + "-centers.txt"),
                    destination_dir / f"level-{level:02d}-centers.txt",
                    expected,
                )
            )
            boundary_source = Path(str(prefix) + "-boundaries.geojson")
            with boundary_source.open(encoding="utf-8") as stream:
                feature_count = len(json.load(stream)["features"])
            if feature_count != expected:
                raise RuntimeError(
                    f"OGC orientation level {level}: expected {expected} boundary features, "
                    f"found {feature_count}"
                )
            files.append(
                copy_checked(
                    boundary_source,
                    destination_dir / f"level-{level:02d}-boundaries.geojson",
                    None,
                )
            )
            runs.append(
                {
                    "level": level,
                    "cell_count": expected,
                    "stderr": stderr.strip(),
                    "stdout_completion": expected_marker(stdout),
                }
            )

    manifest = {
        "schema": 1,
        "source": "DGGRID CLI black-box output",
        "dggrid_version": PINNED_DGGRID_VERSION,
        "dggrid_commit_used_to_build_local_executable": PINNED_DGGRID_COMMIT,
        "dggrid_executable_sha256": sha256(executable),
        "implementation_reference": False,
        "grid": spec.name,
        "purpose": "OGC Annex B orientation cross-check; not an OGC ZIRS oracle",
        "orientation": {
            "vertex_0_longitude_degrees": longitude,
            "vertex_0_spherical_authalic_latitude_degrees": latitude,
            "corresponding_wgs84_geodetic_latitude_degrees": 58.397145907431,
            "vertex_0_azimuth_degrees": 0.0,
        },
        "sphere_radius": 1.0,
        "precision_digits_after_decimal": 15,
        "address_type": "HIERNDX",
        "hierarchy": "Z3",
        "levels": runs,
        "files": sorted(files, key=lambda item: item["path"]),
    }
    (destination_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def expected_marker(stdout: str) -> str:
    matches = re.findall(r"\* generated [0-9][0-9,]* cells", stdout)
    return matches[-1] if matches else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dggrid", type=Path, required=True)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "test" / "oracles",
    )
    arguments = parser.parse_args()
    executable = arguments.dggrid.resolve()
    if not executable.is_file():
        parser.error(f"DGGRID executable does not exist: {executable}")
    generate(executable, arguments.output_root.resolve())


if __name__ == "__main__":
    main()
