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
    """Generate a Gromacs index file with four groups.

    Groups generated:
      - System     (all atoms)
      - lnb_layer  (lipid residues)
      - lnb_gas    (gas residues)
      - solute     (ion + water residues)
    
    Logic:
      1. Read GRO file line by line
      2. Parse atom records (lines 3 to second-to-last)
      3. Extract resname from columns 5-10 (0-based: [5:10])
      4. Classify atoms based on resname into groups
      5. Format output as GROMACS .ndx format (15 indices per line)
    """

    lipids_set: Set[str] = {name.strip().upper() for name in lipid_resnames}
    gases_set: Set[str] = {name.strip().upper() for name in gas_resnames}
    waters_set: Set[str] = {name.strip().upper() for name in water_resnames}
    ions_set: Set[str] = {name.strip().upper() for name in ion_resnames}

    groups: Dict[str, List[int]] = {
        "System": [],
        "lnb_layer": [],
        "lnb_gas": [],
        "solute": [],
    }

    with open(gro_path, "r", encoding="utf-8") as gro_file:
        lines = gro_file.readlines()

    # Atom records are between the 3rd line and before the box vectors line
    # Line 1: title
    # Line 2: number of atoms
    # Lines 3 to N: atom records (N = 2 + atom_count)
    # Line N+1: box vectors (may be followed by empty line)
    # 
    # Determine where atom records end by checking for box vector format
    # Box vectors are typically 3 or 9 floating point numbers
    # Atom records have resname in columns 5-10, which is not all numbers
    atom_count_declared = int(lines[1].strip())
    
    # Check if last line is empty
    if not lines[-1].strip():
        # Last line is empty, box vectors should be second-to-last
        atom_lines = lines[2:-2]
    else:
        # Check if last line looks like box vectors (all numbers, no resname)
        last_line_parts = lines[-1].strip().split()
        if len(last_line_parts) >= 3:
            try:
                # Try to parse as numbers (box vector format)
                float(last_line_parts[0])
                float(last_line_parts[1])
                float(last_line_parts[2])
                # If successful and no resname in expected position, it's box vectors
                if len(lines[-1]) <= 10 or not lines[-1][5:10].strip():
                    atom_lines = lines[2:-1]
                else:
                    # Has resname, it's an atom record
                    atom_lines = lines[2:]
            except ValueError:
                # Not all numbers, it's an atom record
                atom_lines = lines[2:]
        else:
            # Not enough parts, treat as atom record
            atom_lines = lines[2:]
    
    # Validate: atom_lines should match declared count
    if len(atom_lines) != atom_count_declared:
        raise ValueError(
            f"Mismatch: declared {atom_count_declared} atoms, "
            f"found {len(atom_lines)} atom record lines"
        )
    for atom_index, line in enumerate(atom_lines, start=1):
        # GRO file format: resnr(5) resname(5) atomname(5) atomnr(5) x(8) y(8) z(8)
        # resname is at positions 5-10 (0-based: [5:10])
        resname = line[5:10].strip().upper()
        
        # System: all atoms
        groups["System"].append(atom_index)

        # Classify atoms based on resname
        if resname in lipids_set:
            groups["lnb_layer"].append(atom_index)
        elif resname in gases_set:
            groups["lnb_gas"].append(atom_index)
        elif resname in ions_set or resname in waters_set:
            groups["solute"].append(atom_index)

    # Write index file
    os.makedirs(os.path.dirname(ndx_path), exist_ok=True)
    with open(ndx_path, "w", encoding="utf-8") as ndx_file:
        # Write groups in specific order: System, lnb_layer, lnb_gas, solute
        group_order = ["System", "lnb_layer", "lnb_gas", "solute"]
        for group_name in group_order:
            indices = groups.get(group_name, [])
            if not indices:
                continue
            ndx_file.write(f"[ {group_name} ]\n")
            for formatted_line in _format_group(indices):
                ndx_file.write(f"{formatted_line}\n")
            ndx_file.write("\n")
