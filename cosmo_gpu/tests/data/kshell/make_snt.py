"""Convert the repo's OXBASH interactions to KSHELL ``.snt`` with KSHELL's own converter.

The isospin -> proton-neutron conversion (phases, normalizations) is done by KSHELL's
``bin/nushell2snt.py``, so it is independent of ``cosmo_gpu.interaction``. This script only
writes the NuShellX ``.sp`` orbit lists and, for Hamiltonians built from several files (FSU9),
multiplies each converted file by its CoSMo ``<scale>`` and its own mass factor for one mass
number, then sums the files line by line.

usage (from cosmo_gpu/): KSHELL_DIR=~/kshell uv run python tests/data/kshell/make_snt.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from cosmo_gpu.modelspace import ModelSpace, load_hamiltonian_spec, load_model_space

HERE = Path(__file__).resolve().parent
OUT = HERE / "snt"

# (Hamiltonian, CoSMo space, mass numbers for per-nucleus files or None for one general file)
TARGETS: list[tuple[str, str, tuple[int, ...] | None]] = [
    ("usdb", "sd", None),
    ("kb3g", "fp", None),
    ("fsu9", "spsdpf", (20, 22, 24, 28, 33)),
]


def write_sp(space: ModelSpace, path: Path) -> None:
    """NuShellX isospin ``.sp`` file, orbits in CoSMo index order, NuShell n starting at 1."""
    orbits = sorted(space.orbits, key=lambda o: o.index)
    lines = [
        f"! {space.name} written by make_snt.py",
        "t",
        f"{space.core_z + space.core_n} {space.core_z}",
        f"{len(orbits)}",
        f"1 {len(orbits)}",
    ]
    lines += [f"{i + 1} {o.n + 1} {o.l} {o.j2}" for i, o in enumerate(orbits)]
    path.write_text("\n".join(lines) + "\n")


def fix_line_count(interaction: Path, tmp: Path) -> Path:
    """Return the file, or a copy whose header count matches its TBME lines.

    ``nushell2snt.py`` reads exactly the declared number of TBMEs, while CoSMo reads them all
    and only warns. ``fp/kb3g.int`` declares 194 and lists 195, so without this the last
    element (p1/2^2 J=1 T=0) would be dropped silently.
    """
    lines = interaction.read_text().splitlines()
    data = [i for i, ln in enumerate(lines) if ln.strip() and ln.lstrip()[0] not in "!#"]
    head = lines[data[0]].split()
    count = int(float(head[0]))
    n_tbme = len(data) - 1
    if count == 0 or abs(count) == n_tbme:
        return interaction
    print(f"warning: {interaction.name} declares {abs(count)} TBMEs but lists {n_tbme}; fixed")
    head[0] = str(n_tbme if count > 0 else -n_tbme)
    lines[data[0]] = "  ".join(head)
    fixed = tmp / f"fixed_{interaction.name}"
    fixed.write_text("\n".join(lines) + "\n")
    return fixed


def nushell2snt(kshell_bin: Path, sp: Path, interaction: Path, out: Path) -> str:
    """Run KSHELL's converter; answer its scaling prompt (asked only without a mass header)."""
    res = subprocess.run(
        [
            sys.executable,
            "-I",
            str(kshell_bin / "nushell2snt.py"),
            str(sp),
            str(interaction),
            str(out),
        ],
        input="0\n",
        capture_output=True,
        text=True,
        check=True,
        cwd=out.parent,
    )
    return res.stdout


class Snt:
    """The parts of a converted ``.snt`` needed to sum files."""

    def __init__(self, path: Path) -> None:
        rows = []
        self.header: list[str] = []
        for raw in path.read_text().splitlines():
            if raw.lstrip().startswith(("!", "#")):
                continue
            tok = raw.split("!")[0].split()
            if tok:
                rows.append(tok)
        n_p, n_n = int(rows[0][0]), int(rows[0][1])
        self.orbit_rows = rows[: 1 + n_p + n_n]
        pos = 1 + n_p + n_n
        n_1b = int(rows[pos][0])
        self.spe = {(int(r[0]), int(r[1])): float(r[2]) for r in rows[pos + 1 : pos + 1 + n_1b]}
        pos += 1 + n_1b
        head = rows[pos]
        self.method = int(head[1])
        self.a0, self.power = (float(head[2]), float(head[3])) if self.method == 1 else (0.0, 0.0)
        self.tbme: dict[tuple[int, ...], float] = {}
        for r in rows[pos + 1 :]:
            key = tuple(int(x) for x in r[:5])
            self.tbme[key] = self.tbme.get(key, 0.0) + float(r[5])

    def mass_factor(self, mass: int) -> float:
        return (mass / self.a0) ** self.power if self.method == 1 else 1.0


def write_summed(path: Path, parts: list[tuple[Snt, float, float]], comment: list[str]) -> None:
    """Write ``sum_f`` of converted files with s.p. scale and TBME factor, method 0."""
    spe: dict[tuple[int, int], float] = {}
    tbme: dict[tuple[int, ...], float] = {}
    for snt, spe_scale, tbme_factor in parts:
        if snt.orbit_rows != parts[0][0].orbit_rows:
            msg = "files to sum have different orbit lists"
            raise ValueError(msg)
        for k, v in snt.spe.items():
            spe[k] = spe.get(k, 0.0) + spe_scale * v
        for k, v in snt.tbme.items():
            tbme[k] = tbme.get(k, 0.0) + tbme_factor * v
    out = [f"! {c}" for c in comment]
    out += ["  " + "  ".join(r) for r in parts[0][0].orbit_rows]
    out.append(f" {len(spe)} 0")
    out += [f"{i:3d} {j:3d} {v:.17g}" for (i, j), v in sorted(spe.items())]
    out.append(f" {len(tbme)} 0")
    out += [" ".join(f"{x:3d}" for x in k) + f" {v:.17g}" for k, v in tbme.items()]
    path.write_text("\n".join(out) + "\n")


def main() -> int:
    kshell_dir = Path(os.environ.get("KSHELL_DIR", Path.home() / "kshell"))
    kshell_bin = kshell_dir / "bin"
    if not (kshell_bin / "nushell2snt.py").exists():
        print(
            f"KSHELL not found in {kshell_dir}; run tools/kshell/build_kshell.sh", file=sys.stderr
        )
        return 1
    OUT.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)
        for ham, space_name, masses in TARGETS:
            space = load_model_space(space_name)
            spec = load_hamiltonian_spec(ham, space_name)
            sp = tmp / f"{space_name}.sp"
            write_sp(space, sp)
            converted = []
            for k, f in enumerate(spec.files):
                out = tmp / f"{ham}_{k}.snt"
                nushell2snt(kshell_bin, sp, fix_line_count(f.path, tmp), out)
                converted.append((Snt(out), f))
                print(f"{ham}: converted {f.path.name} ({len(converted[-1][0].tbme)} TBMEs)")
            comment = [
                f"{ham} ({space_name}) converted from CoSMo interaction files by KSHELL's",
                "nushell2snt.py (GaffaSnobb/kshell 41ea03a); see tests/data/kshell/make_snt.py",
                "files: "
                + ", ".join(
                    f"{f.path.name} x{f.spe_scale:g}/{f.tbme_scale:g}" for _, f in converted
                ),
            ]
            if masses is None:
                if (
                    len(spec.files) != 1
                    or spec.files[0].spe_scale != 1
                    or spec.files[0].tbme_scale != 1
                ):
                    msg = f"{ham}: a general file needs exactly one unscaled interaction file"
                    raise ValueError(msg)
                text = out.read_text()
                (OUT / f"{ham}_from_int.snt").write_text(
                    "\n".join(f"! {c}" for c in comment) + "\n" + text
                )
                continue
            for mass in masses:
                parts = [(s, f.spe_scale, f.tbme_scale * s.mass_factor(mass)) for s, f in converted]
                write_summed(
                    OUT / f"{ham}_A{mass}.snt", parts, [*comment, f"scaled for A = {mass}"]
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
