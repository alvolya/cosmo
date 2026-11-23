"""Run KSHELL on the validation cases and write ``kshell_reference.json``.

Partitions come from KSHELL's own ``gen_partition`` (imported, not driven through its prompts),
so the Nħω truncation is KSHELL's. Eigenvalues are read at full precision from the ``.wav``
header, which KSHELL writes as a raw stream: ``neig`` (int32), ``mtotal`` (int32), ``neig``
eigenvalues (float64), ``neig`` values of ``jj`` (int32). The log prints only 6 decimals.

usage (from cosmo_gpu/, ideally in a Slurm job, see kshell.sb):
    KSHELL_DIR=~/kshell uv run python tests/data/kshell/run_kshell.py [case-name ...]
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SNT = HERE / "snt"
RUN = HERE / "runs"
REFERENCE = HERE / "kshell_reference.json"
KSHELL_COMMIT = "41ea03a6707fb588e176e906edea85be2a46ac79"
TOL = 1e-10


@dataclass(frozen=True)
class Case:
    """One KSHELL run: nucleus, interaction file, valence numbers, parity and truncation."""

    name: str
    nucleus: str
    snt: str
    z_val: int
    n_val: int
    parity: int
    hw: tuple[int, int] | None = None
    n_eigen: int = 10
    slow: bool = False
    transit: bool = False


CASES = [
    Case("20Ne", "20Ne", "usdb_from_int.snt", 2, 2, 1, transit=True),
    Case("22Na", "22Na", "usdb_from_int.snt", 3, 3, 1),
    Case("24Mg", "24Mg", "usdb_from_int.snt", 4, 4, 1, transit=True),
    Case("25Mg", "25Mg", "usdb_from_int.snt", 4, 5, 1),
    Case("26Al", "26Al", "usdb_from_int.snt", 5, 5, 1, slow=True),
    Case("28Si", "28Si", "usdb_from_int.snt", 6, 6, 1, slow=True),
    Case("44Ti", "44Ti", "kb3g_from_int.snt", 2, 2, 1, transit=True),
    Case("46Ti", "46Ti", "kb3g_from_int.snt", 2, 4, 1, slow=True),
    Case("33Al_0hw+", "33Al", "fsu9_A33.snt", 13, 20, 1, (0, 0), transit=True),
    Case("33Al_1hw-", "33Al", "fsu9_A33.snt", 13, 20, -1, (1, 1)),
    Case("20Ne_1hw-", "20Ne", "fsu9_A20.snt", 10, 10, -1, (1, 1)),
    Case("22Na_1hw-", "22Na", "fsu9_A22.snt", 11, 11, -1, (1, 1), slow=True),
]


def write_partition(kshell_bin: Path, snt: Path, case: Case, out: Path) -> None:
    """Write the ``.ptn`` with KSHELL's ``gen_partition.ModelSpace`` (the steps of its main)."""
    sys.path.insert(0, str(kshell_bin))
    import gen_partition as gp

    with snt.open() as fp:
        n_jorb = [0, 0]
        n_jorb[0], n_jorb[1], _, _ = gp.read_comment_skip(fp)
        norb, lorb, jorb, itorb = [], [], [], []
        for _ in range(sum(n_jorb)):
            arr = gp.read_comment_skip(fp)
            norb.append(arr[1])
            lorb.append(arr[2])
            jorb.append(arr[3])
            itorb.append(arr[4])
    ms = gp.ModelSpace((case.z_val, case.n_val), norb, lorb, jorb, itorb)
    with out.open("w") as fpout:
        if case.hw is not None:
            ms.set_hw_truncation(list(case.hw))
            fpout.write(f"# hw trucnation,  min hw = {case.hw[0]} ,   max hw = {case.hw[1]}\n")
        ms.gen_ptn_pn()
        ms.ptn_combined(case.parity)
        ms.strip_ptn_pn()
        ms.write_ptn_pn(fpout, case.parity, str(snt))
        ms.write_ptn_combined(fpout)


def write_input(case: Case, path: Path, ptn: Path, wav: Path) -> None:
    """KSHELL namelist: thick-restart Lanczos, tight tolerance, no Lawson term of its own."""
    mtot = (case.z_val + case.n_val) % 2
    path.write_text(
        f"""&input
  beta_cm = 0.0
  fn_int = "{SNT / case.snt}"
  fn_ptn = "{ptn.name}"
  fn_save_wave = "{wav.name}"
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 2000
  mode_lv_hdd = 0
  mtot = {mtot}
  n_eigen = {case.n_eigen}
  n_restart_vec = {max(15, case.n_eigen + 5)}
  tol = {TOL}
&end
"""
    )


def read_wav_header(wav: Path) -> tuple[int, list[float], list[int]]:
    """Return (2M, eigenvalues, jj) from the start of a KSHELL ``.wav`` stream file."""
    with wav.open("rb") as f:
        neig, mtotal = np.fromfile(f, dtype="<i4", count=2)
        evals = np.fromfile(f, dtype="<f8", count=int(neig))
        jj = np.fromfile(f, dtype="<i4", count=int(neig))
    return int(mtotal), evals.tolist(), jj.tolist()


#: Transition operators for transit.exe (free-nucleon g factors as CoSMo's XSHLEMB defaults).
EFF_CHARGE = (1.5, 0.5)
GL = (1.0, 0.0)
GS = (2 * 2.79285, 2 * -1.91304)
TRANSIT_ROW = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)"
    r"\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)"
)


def run_transit(kshell_bin: Path, case: Case, work: Path) -> dict[str, object]:
    """Run transit.exe on the case's own wavefunctions; return b, E2 and M1 tables."""
    inp = work / f"{case.name}_tr.input"
    log = work / f"{case.name}_tr.log"
    inp.write_text(
        f"""&input
  fn_int   = "{SNT / case.snt}"
  fn_ptn_l = "{case.name}.ptn"
  fn_ptn_r = "{case.name}.ptn"
  fn_load_wave_l = "{case.name}.wav"
  fn_load_wave_r = "{case.name}.wav"
  hw_type = 1
  eff_charge = {EFF_CHARGE[0]}, {EFF_CHARGE[1]}
  gl = {GL[0]}, {GL[1]}
  gs = {GS[0]}, {GS[1]}
&end
"""
    )
    with log.open("w") as out:
        subprocess.run(
            [str(kshell_bin / "transit.exe"), inp.name],
            cwd=work,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
        )
    text = log.read_text()
    b = re.search(r"b =\s+([\d.]+) fm", text)
    tables: dict[str, list[list[float]]] = {"E2": [], "M1": []}
    current = None
    for line in text.splitlines():
        if line.strip().startswith(("E2 transition", "M1 transition")):
            current = line.split()[0]
            continue
        m = TRANSIT_ROW.match(line)
        if current and m:
            j2f, idf, _, j2i, idi, _, _, mred, b_fwd, _ = m.groups()
            tables[current].append(
                [int(j2f), int(idf), int(j2i), int(idi), float(mred), float(b_fwd)]
            )
        elif current and not line.strip():
            continue
    return {
        "b_fm": float(b.group(1)) if b else None,
        "eff_charge": EFF_CHARGE,
        "gl": GL,
        "gs": GS,
        # rows: [2J_first, idx_first, 2J_second, idx_second, Mred, B(first -> second)]
        "E2": tables["E2"],
        "M1": tables["M1"],
    }


def read_occupations(text: str) -> list[dict[str, list[float]]]:
    """Per state (energy order), the ``<p Nj>`` and ``<n Nj>`` lines of a kshell.exe log."""
    out: list[dict[str, list[float]]] = []
    lines = text.splitlines()
    for k, line in enumerate(lines):
        if "<H>:" in line and "<JJ>:" in line:
            occ: dict[str, list[float]] = {}
            for nxt in lines[k + 1 : k + 6]:
                if nxt.strip().startswith("<p Nj>"):
                    occ["p"] = [float(x) for x in nxt.split()[2:]]
                if nxt.strip().startswith("<n Nj>"):
                    occ["n"] = [float(x) for x in nxt.split()[2:]]
            out.append(occ)
    return out


def run_case(kshell_bin: Path, case: Case) -> dict[str, object]:
    """Run one case in ``runs/<name>/`` and return its reference record."""
    work = RUN / case.name
    work.mkdir(parents=True, exist_ok=True)
    ptn, wav, inp, log = (work / f"{case.name}.{ext}" for ext in ("ptn", "wav", "input", "log"))
    for f in (ptn, wav):
        f.unlink(missing_ok=True)
    write_partition(kshell_bin, SNT / case.snt, case, ptn)
    write_input(case, inp, ptn, wav)
    t0 = time.perf_counter()
    with log.open("w") as out:
        subprocess.run(
            [str(kshell_bin / "kshell.exe"), inp.name],
            cwd=work,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
            env={**os.environ, "GFORTRAN_UNBUFFERED_PRECONNECTED": "y"},
        )
    seconds = time.perf_counter() - t0
    for tmp in work.glob("tmp_*"):
        tmp.unlink()
    text = log.read_text()
    dim = re.search(r"total m-scheme dimension\s+(\d+)", text)
    converged = "H converged" in text
    m2, evals, jj = read_wav_header(wav)
    extra: dict[str, object] = {"occupations": read_occupations(text)}
    if case.transit:
        extra["transit"] = run_transit(kshell_bin, case, work)
    return {
        **asdict(case),
        "dimension": int(dim.group(1)) if dim else None,
        "m2": m2,
        "energies": evals,
        "jj": jj,
        "converged": converged,
        "tol": TOL,
        "seconds": round(seconds, 2),
        "kshell_commit": KSHELL_COMMIT,
        **extra,
    }


def main(names: list[str]) -> int:
    kshell_bin = Path(os.environ.get("KSHELL_DIR", Path.home() / "kshell")) / "bin"
    if not (kshell_bin / "kshell.exe").exists():
        print(
            f"kshell.exe not found in {kshell_bin}; run tools/kshell/build_kshell.sh",
            file=sys.stderr,
        )
        return 1
    records = json.loads(REFERENCE.read_text()) if REFERENCE.exists() else {}
    for case in CASES:
        if names and case.name not in names:
            continue
        rec = run_case(kshell_bin, case)
        records[case.name] = rec
        print(
            f"{case.name}: D={rec['dimension']}  E0={rec['energies'][0]:.10f}  {rec['seconds']} s"
        )
        REFERENCE.write_text(json.dumps(records, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
