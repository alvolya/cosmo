"""Extract CoSMo DEXHVJTcsb energies, spins and convergence into JSON test references.

Run from this directory: ``python extract_cosmo.py``. Writes ``cosmo_usdb_reference.json`` and
``cosmo_fsu9_reference.json``.

CoSMo prints ``State: k E=<energy> JJ=<int(J(J+1))>``; the truncated J(J+1) still fixes J once
the parity of 2J (that of A) is known. After ``writing output to ...HH.EE`` it prints one row per
state: energy, eigenvalue error, eigenvector error. A large eigenvector error marks a state the
Davidson solver did not converge. ``cosmo_fsu9/dex-22Na_1hw-.log`` is the stage-0 benchmark run
(``bench-psdpf/runs/22Na_bin-O3_t16/run.log``, 400 Davidson iterations).
"""

import json
import re
from pathlib import Path

HERE = Path(__file__).parent
NUM = re.compile(r"^[-+.\de]+ [-+.\de]+ [-+.\de]+$")


def parse(log: Path, mass: int) -> dict[str, object]:
    text = log.read_text()
    m = re.search(r"dimension (\d+)", text)
    if m is None:
        raise ValueError(f"{log}: no dimension")
    E, J2 = [], []
    for e, jj in re.findall(r"State: \d+ E= (\S+) JJ= (\S+)", text):
        E.append(float(e))
        j2 = next(
            j2
            for j2 in range(mass % 2, 80, 2)
            if int((j2 / 2) * (j2 / 2 + 1) + 1e-9) == round(float(jj))
        )
        J2.append(j2)
    # first run of len(E) consecutive numeric rows after "writing output to"
    vec_err = [0.0] * len(E)  # dense diagonalization prints no table
    tail = text.split("writing output to", 1)[-1].splitlines()
    run: list[str] = []
    for line in tail:
        if NUM.match(line.strip()):
            run.append(line)
            if len(run) == len(E):
                vec_err = [float(r.split()[2]) for r in run]
                break
        else:
            run = []
    return {"dim": int(m.group(1)), "energies": E, "j2": J2, "vector_error": vec_err}


def split_combined(out: Path, dest: Path) -> None:
    """Split the combined first reference job output into per-nucleus logs."""
    for chunk in out.read_text().split("== ")[1:]:
        nuc = chunk.split("_")[0]
        if nuc in ("24Mg", "25Mg", "28Si"):
            (dest / f"dex-{nuc}.log").write_text(chunk)


def main() -> None:
    usdb = HERE / "cosmo_usdb"
    split_combined(usdb / "ref-21582054.out", usdb)
    refs = {
        "cosmo_usdb_reference.json": {
            nuc: (usdb / f"dex-{nuc}.log", int(re.match(r"\d+", nuc).group(0)))  # type: ignore[union-attr]
            for nuc in ("18O", "19O", "20Ne", "22Na", "24Mg", "25Mg", "26Al", "28Si")
        },
        "cosmo_fsu9_reference.json": {
            key: (HERE / "cosmo_fsu9" / f"dex-{key}.log", int(key[:2]))
            for key in ("33Al_0hw+", "33Al_1hw-", "20Ne_1hw-", "22Na_1hw-")
        },
    }
    for name, systems in refs.items():
        out = {}
        for key, (log, mass) in systems.items():
            out[key] = parse(log, mass)
            r = out[key]
            print(name, key, r["dim"], r["energies"][:2], r["j2"][:4], max(r["vector_error"]))
        (HERE / name).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
