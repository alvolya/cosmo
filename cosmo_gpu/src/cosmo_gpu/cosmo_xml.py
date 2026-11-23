"""CoSMo system XML files: the nucleus, space, Hamiltonian and truncations of one CoSMo job.

A system file (as written by ``cosmoxml`` or the cosmo skill's ``make_system.py``) holds an
embedded ``<valencespace>`` and ``<hamiltonian>`` plus one or more ``<system>`` elements::

    <system Zv="10" Nv="10" Jz="0" P="-" name="20Ne_fsu9_M0T0-">
        <hwrejection HWMAX="1" HWMIN="1"/>
        <rejection NMIN="0" NMAX="2"> <orbital index="3" type="pn"/> ... </rejection>
        <?ManyBodyStates name="....mbs" size="26600" creator="sysmbsOMP"?>
        <state J="2" P="-" E="-179.396" name="2-(1)"> ... </state>
    </system>

Rejection orbitals are matched by their CoSMo ``index``, and ``type`` selects protons (``p``),
neutrons (``n``) or both (``pn``), as in ``XsysmbsOMP``. CoSMo's ``<matrix JJ=...>`` J^2 shift
is not supported.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from cosmo_gpu.modelspace import (
    HamiltonianSpec,
    ModelSpace,
    parse_hamiltonian,
    parse_valence_space,
    spin2,
)
from cosmo_gpu.truncation import OccupationRule


@dataclass(frozen=True)
class CosmoState:
    """A level stored by CoSMo in the system file (energies to 3 decimals)."""

    name: str
    j2: int
    parity: int
    energy: float


@dataclass(frozen=True)
class CosmoSystem:
    """One CoSMo system: everything needed to rebuild its basis and Hamiltonian."""

    name: str
    space: ModelSpace
    hamiltonian: HamiltonianSpec
    z_val: int
    n_val: int
    m2: int
    parity: int
    hw: tuple[int, int] | None
    occupation: tuple[OccupationRule, ...]
    cosmo_dimension: int | None
    cosmo_states: tuple[CosmoState, ...]


def _parse(path: Path) -> ET.Element:
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_pis=True))
    return ET.parse(path, parser).getroot()


def _parity(text: str) -> int:
    if text not in ("+", "-"):
        msg = f"bad parity {text!r}"
        raise ValueError(msg)
    return 1 if text == "+" else -1


def _rejection(rej: ET.Element, space: ModelSpace) -> OccupationRule:
    proton = [False] * space.n_orbits
    neutron = [False] * space.n_orbits
    for o in rej.findall("orbital"):
        pos = space.orbit_position(int(o.get("index", "0")))
        kind = o.get("type", "pn")
        if kind not in ("p", "n", "pn"):
            msg = f"bad rejection orbital type {kind!r}"
            raise ValueError(msg)
        if "p" in kind:
            proton[pos] = True
        if "n" in kind:
            neutron[pos] = True
    return OccupationRule(
        proton_orbits=tuple(proton),
        neutron_orbits=tuple(neutron),
        nmin=int(rej.get("NMIN", "0")),
        nmax=int(rej.get("NMAX", "0")),
    )


def _dimension(sysnode: ET.Element) -> int | None:
    for child in sysnode:
        tag: object = child.tag  # processing instructions carry the factory as their tag
        if tag is ET.ProcessingInstruction and (child.text or "").startswith("ManyBodyStates"):
            m = re.search(r'size="(\d+)"', child.text or "")
            if m:
                return int(m.group(1))
    return None


def system_names(path: Path) -> list[str]:
    """Names of the ``<system>`` elements in a CoSMo XML file."""
    return [s.get("name", "") for s in _parse(path).findall("system")]


def read_cosmo_xml(
    path: Path, system: str | None = None, interactions_dir: Path | None = None
) -> CosmoSystem:
    """Read one system from a CoSMo XML file.

    Parameters
    ----------
    path
        The system file.
    system
        Name of the ``<system>`` to read; may be omitted when the file holds only one.
    interactions_dir
        Where interaction files are searched (default: as ``find_interaction_file``).

    Raises
    ------
    ValueError
        If the file lacks a valence space, Hamiltonian or the requested system.
    NotImplementedError
        For CoSMo's ``<matrix JJ=...>`` J^2 shift.
    """
    root = _parse(path)
    vs = root.find("valencespace")
    ham = root.find("hamiltonian")
    if vs is None or ham is None:
        msg = f"{path.name}: needs an embedded <valencespace> and <hamiltonian>"
        raise ValueError(msg)
    systems = root.findall("system")
    if system is not None:
        systems = [s for s in systems if s.get("name") == system]
    if len(systems) != 1:
        names = [s.get("name") for s in root.findall("system")]
        msg = f"{path.name}: choose one system (found {names})"
        raise ValueError(msg)
    sysnode = systems[0]
    for m in sysnode.findall("matrix"):
        if float(m.get("JJ", "0")) != 0.0:
            msg = f"{path.name}: the <matrix JJ=...> J^2 shift is not supported"
            raise NotImplementedError(msg)
    space = parse_valence_space(vs, interactions_dir)
    hw_node = sysnode.find("hwrejection")
    hw = (
        (int(hw_node.get("HWMIN", "0")), int(hw_node.get("HWMAX", "0")))
        if hw_node is not None
        else None
    )
    states = tuple(
        CosmoState(
            name=st.get("name", ""),
            j2=spin2(st.get("J", "0")),
            parity=_parity(st.get("P", "+")),
            energy=float(st.get("E", "nan")),
        )
        for st in sysnode.findall("state")
    )
    return CosmoSystem(
        name=sysnode.get("name", ""),
        space=space,
        hamiltonian=parse_hamiltonian(ham, interactions_dir),
        z_val=int(sysnode.get("Zv", "0")),
        n_val=int(sysnode.get("Nv", "0")),
        m2=spin2(sysnode.get("Jz", "0")),
        parity=_parity(sysnode.get("P", "+")),
        hw=hw,
        occupation=tuple(_rejection(r, space) for r in sysnode.findall("rejection")),
        cosmo_dimension=_dimension(sysnode),
        cosmo_states=states,
    )
