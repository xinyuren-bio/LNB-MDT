import os
from typing import Iterable, List, Dict, Set


def _format_group(indices: Iterable[int]) -> List[str]:
    """Format atom indices into Gromacs .ndx style lines (15 integers per line)."""
    buffer: List[str] = []
    lines: List[str] = []
    for idx in indices:
        buffer.append(str(idx))
        if len(buffer) == 15:
            lines.append(" ".join(buffer))
            buffer = []
    if buffer:
        lines.append(" ".join(buffer))
    return lines


def generate_ndx(
    gro_path: str,
    ndx_path: str,
    lipid_resnames: Iterable[str],
    gas_resnames: Iterable[str],
    water_resnames: Iterable[str],
    ion_resnames: Iterable[str],
) -> None:
    """Generate a basic Gromacs index file with standard groups.

    Groups generated:
      - System (all atoms)
      - IONS   (ion residues)
      - Lipids (lipid residues)
      - Gas    (gas residues)
      - LNB    (lipids + gas residues)
      - W      (water residues)
    """

    lipids_set: Set[str] = {name.strip().upper() for name in lipid_resnames}
    gases_set: Set[str] = {name.strip().upper() for name in gas_resnames}
    waters_set: Set[str] = {name.strip().upper() for name in water_resnames}
    ions_set: Set[str] = {name.strip().upper() for name in ion_resnames}

    groups: Dict[str, List[int]] = {
        "System": [],
        "IONS": [],
        "Lipids": [],
        "Gas": [],
        "W": [],
    }

    with open(gro_path, "r", encoding="utf-8") as gro_file:
        lines = gro_file.readlines()

    # Atom records are between the 3rd line and the last line (box vectors)
    atom_lines = lines[2:-1]
    for atom_index, line in enumerate(atom_lines, start=1):
        resname = line[5:10].strip().upper()
        groups["System"].append(atom_index)

        if resname in ions_set:
            groups["IONS"].append(atom_index)
        if resname in lipids_set:
            groups["Lipids"].append(atom_index)
        if resname in gases_set:
            groups["Gas"].append(atom_index)
        if resname in waters_set:
            groups["W"].append(atom_index)

    # LNB group combines lipids and gas residues
    lnb_indices = sorted(set(groups["Lipids"]) | set(groups["Gas"]))
    if lnb_indices:
        groups["LNB"] = lnb_indices

    os.makedirs(os.path.dirname(ndx_path), exist_ok=True)
    with open(ndx_path, "w", encoding="utf-8") as ndx_file:
        for group_name, indices in groups.items():
            if not indices:
                continue
            ndx_file.write(f"[ {group_name} ]\n")
            for formatted_line in _format_group(indices):
                ndx_file.write(f"{formatted_line}\n")
            ndx_file.write("\n")
