# Example Lipid Nanobubble (LNB) System

This directory contains example molecular dynamics simulation files for a lipid nanobubble system.

## System Description

### Topology File (`lnb.gro`)

The topology file `lnb.gro` represents a lipid nanobubble system simulated using the **Martini 3.0** coarse-grained force field. The system composition is as follows:

- **Lipid composition**: DPPC (1,2-dipalmitoyl-sn-glycero-3-phosphocholine) : DAPC (1,2-diarachidoyl-sn-glycero-3-phosphocholine) : CHOL (cholesterol) = **5 : 3 : 2** (molar ratio)
- **Gas molecules**: O₂ (oxygen) molecules are included within the nanobubble
- **Solvent**: Water molecules and ions have been removed from the topology file to reduce file size

### Trajectory File (`lnb.xtc`)

The trajectory file `lnb.xtc` contains the molecular dynamics trajectory data extracted from the simulation time period of **50-60 ns**. This time window was selected to provide a representative sample of the system dynamics while maintaining a manageable file size for distribution and analysis purposes.

## Notes

- The system uses the Martini 3.0 coarse-grained force field
- Solvent and ions were removed post-simulation to optimize file storage
- The trajectory represents a 10 ns time window from the production run

