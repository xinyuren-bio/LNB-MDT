#!/usr/bin/env python3
"""
LNB System Generator with Automatic Folder Structure
====================================================

This script generates a lipid nanobubble system and automatically creates
a complete folder structure with all necessary files for GROMACS simulations.

Usage:
    python preparation/lnb_generate_with_folder.py [OPTIONS] -o <folder_name>

The -o parameter specifies the output folder name (not file path).
Generated files will be automatically named (system.gro, system.top, system.ndx).

Example:
    python preparation/lnb_generate_with_folder.py \
        -d 1 -r 5 -x 20 -y 20 -z 20 \
        -u DPPC:50 -u DAPC:30 -u CHOL:20 \
        -sol W -salt 0.15 -a 1 \
        -gas O2 -gden 200 \
        -o my_lnb_system

This will create a folder 'LNB_MDT_my_lnb_system' containing:
    - system.gro
    - system.top
    - system.ndx
    - parameter.txt
    - toppar/ (with all .itp files)
    - *.mdp files (simulation parameters)
"""

import os
import sys
import subprocess
import shutil
from datetime import datetime
from pathlib import Path

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LNB_GENER_SCRIPT = os.path.join(SCRIPT_DIR, 'lnb_gener_martini3.py')
FILES_DIR = os.path.join(SCRIPT_DIR, 'files')


def parse_arguments():
    """Parse command line arguments, similar to lnb_gener_martini3.py"""
    args = sys.argv[1:]
    
    # Check for help
    if '-h' in args or '--help' in args:
        print(__doc__)
        print("\nAll parameters are the same as lnb_gener_martini3.py, except:")
        print("  -o: Output folder name (not file path)")
        print("\nFor detailed parameter descriptions, see:")
        print("  python preparation/lnb_gener_martini3.py -h")
        sys.exit(0)
    
    # Auto-fix common mistakes: split concatenated values (e.g., "20-u" -> "20", "-u")
    fixed_args = []
    i = 0
    while i < len(args):
        arg = args[i]
        # Check if this is a value (not starting with -) that might be concatenated
        if not arg.startswith('-') and i > 0:
            prev_arg = args[i - 1]
            # Check if value contains a dash followed by a letter (likely a missing space)
            if '-' in arg and len(arg) > 1:
                # Try to split: look for pattern like "20-u" or "50-DPPC"
                parts = arg.split('-', 1)
                if len(parts) == 2 and parts[0] and parts[1]:
                    # Check if first part looks like a number
                    try:
                        float(parts[0])
                        # First part is a number, second part likely starts a new parameter
                        fixed_args.append(parts[0])
                        fixed_args.append('-' + parts[1])
                        i += 1
                        continue
                    except ValueError:
                        pass
        fixed_args.append(arg)
        i += 1
    
    args = fixed_args
    
    # Find -o parameter and extract folder name
    output_folder = None
    if '-o' in args:
        idx = args.index('-o')
        if idx + 1 < len(args):
            output_folder = args[idx + 1]
            # Remove -o and its value from args for passing to lnb_gener_martini3.py
            args = args[:idx] + args[idx+2:]
        else:
            print("Error: -o requires a folder name")
            sys.exit(1)
    else:
        print("Error: -o parameter (output folder name) is required")
        print("Usage: python preparation/lnb_generate_with_folder.py [OPTIONS] -o <folder_name>")
        sys.exit(1)
    
    return args, output_folder


def create_temporary_gro_path(output_folder):
    """Create a temporary GRO file path for lnb_gener_martini3.py"""
    # Use a temporary name that will be renamed later
    # Use absolute path to ensure we can find the generated ndx file
    cwd = os.getcwd()
    temp_gro = os.path.join(cwd, f"{output_folder}_temp.gro")
    return temp_gro


def run_lnb_generator(args, temp_gro_path):
    """Run lnb_gener_martini3.py with modified arguments"""
    # Build command
    command = [sys.executable, LNB_GENER_SCRIPT] + args + ['-o', temp_gro_path]
    
    # Run the generator
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running lnb_gener_martini3.py:")
        print(result.stderr)
        print(f"\nCommand that failed:")
        print(f"  {' '.join(command)}")
        print(f"\nTip: Make sure all parameters have proper spacing.")
        print(f"     For example: -z 20 -u DPPC:50 (not -z 20-u DPPC:50)")
        sys.exit(1)
    
    # Extract molecule information from stderr (standard output format)
    molecule_info = result.stderr.strip() if result.stderr else ""
    
    return molecule_info


def create_topology_file(output_dir, molecule_info):
    """Create system.top file with proper includes"""
    top_content = """;
;  
; Example topology file for MARTINI 3 
;  

; First include the file containing all particle definitions,  
; the interaction matrix, plus the topology for water.  


; Then include the file(s) containing the topologies of other  
; molecules present in your system.  

#include "toppar/martini_v3.0.0.itp"
#include "toppar/martini_v3.0.0_ffbonded.itp"
#include "toppar/martini_v3.0.0_ions.itp"  
#include "toppar/martini_v3.0.0_phospholipids.itp"
#include "toppar/martini_v3.0.0_gas.itp" 
#include "toppar/martini_v3.0.0_solvents.itp"


; Define a name for your system  

[ system ]  
Lipid-Nanobubble  

; Define the composition of your system  
; The molecule names should correspond to those defined in the itp file(s).  

[ molecules ]
{}  

""".format(molecule_info)
    
    top_path = os.path.join(output_dir, 'topol.top')
    with open(top_path, 'w') as f:
        f.write(top_content)
    
    return top_path


def create_parameter_file(main_folder_path, args_dict, molecule_info):
    """Create parameter.txt file with all generation parameters"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    parameter_content = f"""LNB-MDT Generation Parameters
===============================

Generation Time: {timestamp}

System Parameters:
-----------------
Box Dimensions:
  - X: {args_dict.get('x', 'N/A')}
  - Y: {args_dict.get('y', 'N/A')}
  - Z: {args_dict.get('z', 'N/A')}
  - D: {args_dict.get('d', '1')}

LNB Parameters:
---------------
Gas Type: {args_dict.get('gas', 'N/A')}
Gas Density: {args_dict.get('gden', 'N/A')}
Area: {args_dict.get('a', 'N/A')}
Radius: {args_dict.get('r', 'N/A')}

Solvent Parameters:
------------------
Solvent Type: {args_dict.get('sol', 'N/A')}
Salt Concentration: {args_dict.get('salt', 'N/A')}

Lipid Composition:
-----------------"""
    
    # Extract lipid information from -u parameters
    lipids = args_dict.get('u', [])
    if lipids:
        for lipid in lipids:
            if ':' in lipid:
                name, count = lipid.split(':', 1)
                parameter_content += f"\n  - {name}: {count}"
    else:
        parameter_content += "\n  - No lipids specified"
    
    parameter_content += f"""

Molecule Information:
-------------------
{molecule_info}

Output Information:
------------------
Output Folder: {os.path.basename(main_folder_path)}
Base Filename: system

Generation Command:
------------------"""
    
    # Reconstruct command
    command_parts = ["python", "preparation/lnb_generate_with_folder.py"]
    for key, value in args_dict.items():
        if key in ['u', 'sol', 'gas'] and isinstance(value, list):
            for v in value:
                command_parts.extend([f'-{key}', str(v)])
        elif value and key != 'o':
            if value is not True:  # Skip boolean flags
                command_parts.extend([f'-{key}', str(value)])
            else:
                command_parts.append(f'-{key}')
    command_parts.extend(['-o', os.path.basename(main_folder_path)])
    
    parameter_content += f"\n{' '.join(command_parts)}\n"
    
    parameter_content += """
Notes:
------
- This file contains all parameters used for system generation
- Generated by LNB-MDT v1.0
- For questions or support, please refer to the documentation
"""
    
    parameter_file_path = os.path.join(main_folder_path, "parameter.txt")
    with open(parameter_file_path, 'w', encoding='utf-8') as f:
        f.write(parameter_content)
    
    return parameter_file_path


def create_folder_structure(output_folder, temp_gro_path, args_dict, molecule_info):
    """Create complete folder structure with all files"""
    # Get current working directory
    cwd = os.getcwd()
    
    # Create main folder name
    main_folder_name = f"LNB_MDT_{output_folder}"
    main_folder_path = os.path.join(cwd, main_folder_name)
    
    # Create main folder
    os.makedirs(main_folder_path, exist_ok=True)
    
    # Create toppar subfolder
    toppar_folder_path = os.path.join(main_folder_path, "toppar")
    os.makedirs(toppar_folder_path, exist_ok=True)
    
    try:
        # Move and rename gro file
        if os.path.exists(temp_gro_path):
            system_gro_path = os.path.join(main_folder_path, "system.gro")
            shutil.copy2(temp_gro_path, system_gro_path)
            os.remove(temp_gro_path)
        else:
            print(f"Warning: {temp_gro_path} not found")
        
        # Create and move top file
        create_topology_file(cwd, molecule_info)
        top_path = os.path.join(cwd, "topol.top")
        if os.path.exists(top_path):
            system_top_path = os.path.join(main_folder_path, "system.top")
            shutil.copy2(top_path, system_top_path)
            os.remove(top_path)
        
        # Handle NDX file - try to find it first, then generate if needed
        system_ndx_path = os.path.join(main_folder_path, "system.ndx")
        found_ndx = False
        
        # Try multiple possible paths where ndx might have been generated
        possible_ndx_paths = [
            os.path.join(cwd, f"{os.path.basename(os.path.splitext(temp_gro_path)[0])}.ndx"),
            os.path.join(cwd, f"{output_folder}_temp.ndx"),
            os.path.splitext(temp_gro_path)[0] + ".ndx",
            os.path.join(cwd, os.path.basename(temp_gro_path).replace('.gro', '.ndx')),
        ]
        
        for alt_path in possible_ndx_paths:
            if os.path.exists(alt_path):
                shutil.copy2(alt_path, system_ndx_path)
                os.remove(alt_path)
                found_ndx = True
                print("✓ NDX file found and moved")
                break
        
        # Generate ndx file manually if it wasn't created by lnb_gener_martini3.py
        if not found_ndx:
            print("Generating NDX file manually...")
            try:
                from preparation.ndx_generator import generate_ndx
                import MDAnalysis as mda
                
                # Load gro to get resnames
                gro_file = os.path.join(main_folder_path, "system.gro")
                u = mda.Universe(gro_file)
                all_resnames = set(u.atoms.resnames)
                
                # Classify resnames - use comprehensive lipid list
                try:
                    from preparation.lipidsInfo_martini3 import ALL_P
                    lipid_keywords = set(ALL_P)
                except:
                    # Fallback to common lipids
                    lipid_keywords = set(['DPPC', 'DAPC', 'CHOL', 'POPC', 'DOPC', 'POPE', 'DOPE', 'POPS', 'DOPS', 
                                         'POPG', 'DOPG', 'DLPC', 'DMPC', 'DSPC', 'DLPE', 'DMPE', 'DSPE', 
                                         'PA', 'PS', 'PG', 'PE', 'PC', 'SM'])
                
                gas_keywords = {"CO2", "N2", "O2", "H2", "AIR"}
                water_keywords = {"W", "PW", "SPC", "SPCM", "FG4W", "FG4W-MS", "BMW"}
                ion_keywords = {"NA", "CL", "Mg", "K"}
                
                lipid_resnames = [r for r in all_resnames if r in lipid_keywords]
                gas_resnames = [r for r in all_resnames if r in gas_keywords]
                water_resnames = [r for r in all_resnames if r in water_keywords]
                ion_resnames = [r for r in all_resnames if r in ion_keywords]
                
                generate_ndx(
                    gro_path=gro_file,
                    ndx_path=system_ndx_path,
                    lipid_resnames=lipid_resnames,
                    gas_resnames=gas_resnames,
                    water_resnames=water_resnames,
                    ion_resnames=ion_resnames,
                )
                print("✓ NDX file generated successfully")
            except Exception as e:
                print(f"Warning: Could not generate NDX file: {e}")
                import traceback
                traceback.print_exc()
        
        # Copy files from files/ directory
        if os.path.exists(FILES_DIR):
            for file in os.listdir(FILES_DIR):
                source_path = os.path.join(FILES_DIR, file)
                if os.path.isfile(source_path):
                    if file.endswith('.itp'):
                        # .itp files go to toppar folder
                        dest_path = os.path.join(toppar_folder_path, file)
                        shutil.copy2(source_path, dest_path)
                    else:
                        # Other files (README, .mdp, etc.) go to main folder
                        dest_path = os.path.join(main_folder_path, file)
                        shutil.copy2(source_path, dest_path)
        
        # Create parameter file
        create_parameter_file(main_folder_path, args_dict, molecule_info)
        
        return main_folder_path
        
    except Exception as e:
        print(f"Error creating folder structure: {e}")
        import traceback
        traceback.print_exc()
        return None


def parse_args_to_dict(original_args, output_folder):
    """Parse arguments list into a dictionary for parameter file"""
    args_dict = {}
    # Reconstruct full args list for parsing
    full_args = original_args + ['-o', output_folder]
    
    i = 0
    while i < len(full_args):
        if full_args[i].startswith('-'):
            key = full_args[i][1:]  # Remove leading '-'
            if i + 1 < len(full_args) and not full_args[i + 1].startswith('-'):
                value = full_args[i + 1]
                if key == 'u':
                    if 'u' not in args_dict:
                        args_dict['u'] = []
                    args_dict['u'].append(value)
                elif key == 'sol':
                    if 'sol' not in args_dict:
                        args_dict['sol'] = []
                    args_dict['sol'].append(value)
                elif key == 'gas':
                    if 'gas' not in args_dict:
                        args_dict['gas'] = []
                    args_dict['gas'].append(value)
                else:
                    args_dict[key] = value
                i += 2
            else:
                args_dict[key] = True
                i += 1
        else:
            i += 1
    return args_dict


def main():
    """Main function"""
    # Parse arguments
    args, output_folder = parse_arguments()
    
    # Create temporary GRO file path
    temp_gro_path = create_temporary_gro_path(output_folder)
    
    # Parse args to dictionary for parameter file
    args_dict = parse_args_to_dict(args, output_folder)
    
    # Run LNB generator
    print(f"Generating LNB system...")
    molecule_info = run_lnb_generator(args, temp_gro_path)
    
    # Create folder structure
    print(f"Creating folder structure...")
    main_folder_path = create_folder_structure(output_folder, temp_gro_path, args_dict, molecule_info)
    
    if main_folder_path:
        print(f"\n✓ Success! System generated in folder:")
        print(f"  {main_folder_path}")
        print(f"\nFolder structure:")
        print(f"  - system.gro")
        print(f"  - system.top")
        print(f"  - system.ndx")
        print(f"  - parameter.txt")
        print(f"  - toppar/ (with all .itp files)")
        print(f"  - *.mdp files (simulation parameters)")
    else:
        print("\n✗ Error: Failed to create folder structure")
        sys.exit(1)


if __name__ == "__main__":
    main()

