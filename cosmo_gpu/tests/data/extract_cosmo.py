"""Extract CoSMo DEXHVJTcsb energies and spins into ../cosmo_usdb_reference.json.

CoSMo prints ``State: k E= <energy> JJ= <int(J(J+1))>``; J(J+1) is truncated to an integer,
which still determines J uniquely.
"""

import json
import re
from pathlib import Path

here = Path(__file__).parent
logs = {"18O": "dex-18O.log", "19O": "dex-19O.log", "20Ne": "dex-20Ne.log"}
big = (here / "ref-21582054.out").read_text().split("== ")
for chunk in big[1:]:
    nuc = chunk.split("_")[0]
    if nuc in ("24Mg", "25Mg", "28Si"):
        (here / f"dex-{nuc}.log").write_text(chunk)
        logs[nuc] = f"dex-{nuc}.log"

out = {}
for nuc, log in logs.items():
    text = (here / log).read_text()
    dim = int(re.search(r"dimension (\d+)", text).group(1))
    E, J2 = [], []
    for e, jj in re.findall(r"State: \d+ E= (\S+) JJ= (\S+)", text):
        E.append(float(e))
        jjf = float(jj)
        # smallest 2J with floor(J(J+1)) == JJ
        odd = int(re.match(r"\d+", nuc).group(0)) % 2
        j2 = next(
            j2 for j2 in range(odd, 80, 2) if int((j2 / 2) * (j2 / 2 + 1) + 1e-9) == round(jjf)
        )
        J2.append(j2)
    # convergence table after "states <dim>": energy, eigenvalue error, eigenvector error
    vec_err = [0.0] * len(E)
    num = r"[-+.\de]+"
    table = re.search(rf"^states {dim}\n((?:{num} {num} {num}\n)+)", text, flags=re.M)
    if table:
        vec_err = [float(row.split()[2]) for row in table.group(1).splitlines()]
    out[nuc] = {"dim": dim, "energies": E, "j2": J2, "vector_error": vec_err}
    print(nuc, dim, len(E), E[:3], J2[:5], max(vec_err))
(here.parent / "cosmo_usdb_reference.json").write_text(json.dumps(out, indent=1) + "\n")
