"""Collect CoSMo's stage-2 operator results into ``cosmo_ops_reference.json``.

Inputs (written by ``ops.sb`` and the ``XSHLSF`` run documented in ``README``-style comments of
``ops.sb``): ``<X>_E2.xml`` and ``<X>_M1.xml`` from ``XSHLEMB`` (``B`` in oscillator units for
E2 and in mu_N^2 for M1), ``<X>_occ.txt`` from ``XSHLAO`` (an XML file), ``dex-<X>.log`` for the
Davidson eigenvector errors, and ``sf-25Mg-24Mg.log`` from
``XSHLSF 25Mg_usdb+.xml -f 24Mg_usdb+.xml`` (per-orbit C2S for neutron removal).

usage (from this directory): python3 extract_ops.py
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
NUCLEI = ("20Ne", "22Na", "24Mg", "25Mg")


def eigenvector_errors(log: Path) -> list[float]:
    """``dV`` of each state, in order, from the lines after ``writing output to``."""
    lines = log.read_text().splitlines()
    start = next(i for i, ln in enumerate(lines) if "writing output to" in ln)
    out = []
    for ln in lines[start + 1 :]:
        tok = ln.split()
        if len(tok) != 3:
            break
        try:
            out.append(float(tok[2]))
        except ValueError:
            break
    return out


def transitions(path: Path) -> list[list[object]]:
    """``[initial name, final name, B]`` for every transition in an ``XSHLEMB`` output."""
    rows: list[list[object]] = []
    for st in ET.parse(path).getroot().iter("state"):
        rows.extend(
            [st.get("name"), tr.get("final"), float(tr.get("B", "nan"))]
            for tr in st.findall("transition")
        )
    return rows


def occupations(path: Path) -> dict[str, dict[str, list[float]]]:
    """Per state name, per orbit name: ``[protons, neutrons]``."""
    out: dict[str, dict[str, list[float]]] = {}
    for st in ET.parse(path).getroot().iter("state"):
        out[st.get("name", "")] = {
            o.get("name", ""): [float(o.get("Z", "nan")), float(o.get("N", "nan"))]
            for o in st.findall("occupation")
        }
    return out


def spectroscopic_factors(log: Path) -> list[list[object]]:
    """``[parent, daughter, orbit (n, l, 2j), C2S]`` from an ``XSHLSF`` log."""
    rows: list[list[object]] = []
    parent = daughter = None
    orbit_re = re.compile(r"^(\d)([spdfgh])(\d+)/2 C2S=(\S+)")
    lmap = "spdfgh"
    for ln in log.read_text().splitlines():
        if ln.startswith("Parent:"):
            parent = ln.split()[1]
        elif ln.startswith("Daughter:"):
            daughter = ln.split()[1]
        else:
            m = orbit_re.match(ln.strip())
            if m and parent and daughter:
                n, l, j2, c2s = int(m.group(1)), lmap.index(m.group(2)), int(m.group(3)), m.group(4)
                rows.append([parent, daughter, [n, l, j2], float(c2s)])
    return rows


def main() -> None:
    """Write ``cosmo_ops_reference.json`` next to this script."""
    ref: dict[str, object] = {
        "e2_charges": [1.5, 0.5],
        "m1": {"gl": [1.0, 0.0], "mu": [2.79285, -1.91304]},
    }
    nuclei: dict[str, object] = {}
    for x in NUCLEI:
        root = ET.parse(HERE / f"{x}_E2.xml").getroot()
        states = [
            {"name": st.get("name"), "E": float(st.get("E", "nan"))} for st in root.iter("state")
        ]
        for s, dv in zip(states, eigenvector_errors(HERE / f"dex-{x}.log"), strict=False):
            s["vector_error"] = dv
        nuclei[x] = {
            "states": states,
            "E2": transitions(HERE / f"{x}_E2.xml"),
            "M1": transitions(HERE / f"{x}_M1.xml"),
            "occupations": occupations(HERE / f"{x}_occ.txt"),
        }
    ref["nuclei"] = nuclei
    ref["sf_25Mg_24Mg"] = spectroscopic_factors(HERE / "sf-25Mg-24Mg.log")
    (HERE.parent / "cosmo_ops_reference.json").write_text(json.dumps(ref, indent=1) + "\n")


if __name__ == "__main__":
    main()
