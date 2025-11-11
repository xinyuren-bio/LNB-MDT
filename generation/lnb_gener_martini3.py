import sys,math,random,cmath,copy
from ndx_generator import generate_ndx

version  = "20211229.LNB"

# This script was written based on the open source script insane.py (Djurre H. de Jong et al. J. Chem. Theory Comput. 2013, 9, 687-697.) by Yuan He and Yuxuan Wang in Dr. Xubo Lin's group@Beihang University (https://linxubo.github.io).
# This script can be used to generate lipid nanobubble (spherical lipid monolayer) using one-line command such as "python3 LNB.py -d 1 -r 5 -x 20 -y 20 -z 20 -u DPPC:50 -u DPPE:50 -sol W -salt 0.15 -a 1 -gas N2 -gden 200 -o output.gro".
# Details about the parameters for the one-line command can be found using "python3 LNB.py -h".
# This script can only be used with python 3

# PROTOLIPID (diacylglycerol), 18 beads
#
# 1-3-4-6-7--9-10-11-12-13-14
#  \| |/  |
#   2 5   8-15-16-17-18-19-20
#

lipidsx = {}
lipidsy = {}
lipidsz = {}
lipidsa = {}
#
# Diacyl glycerols
moltype = "lipids"
lipidsx[moltype] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipidsy[moltype] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipidsa.update({

    "DAPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  D1B  D2B  D3B  C4B    -    -   "),
    "DAPC": (moltype, "    -    -  NC3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  C4A    -    -  D1B  D2B  D3B  D4B  C4B    - "),
    "DAPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  D1B  D2B  D3B  C4B    -    -   "),
    "DAPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  D1B  D2B  D3B  C4B    -    -   "),
    "DAPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  D1B  D2B  D3B  C4B    -    -   "),
    "DBPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DBPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DBPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DBPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DBPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "DFPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  D3A  D4A    -    -    -  C1B  D2B  D3B  D4B    -    -   "),
    "DFPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  D3A  D4A    -    -    -  C1B  D2B  D3B  D4B    -    -   "),
    "DFPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  D3A  D4A    -    -    -  C1B  D2B  D3B  D4B    -    -   "),
    "DFPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  D3A  D4A    -    -    -  C1B  D2B  D3B  D4B    -    -   "),
    "DFPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  D3A  D4A    -    -    -  C1B  D2B  D3B  D4B    -    -   "),

    "DGPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DGPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DGPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DGPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DGPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),

    "DIPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  D2B  D3B  C4B    -    -   "),
    "DIPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  D2B  D3B  C4B    -    -   "),
    "DIPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  D2B  D3B  C4B    -    -   "),
    "DIPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  D2B  D3B  C4B    -    -   "),
    "DIPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  D2B  D3B  C4B    -    -   "),
    "DLPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A    -    -    -    -  C1B  C2B  C3B    -    -    -   "),

    "DLPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A    -    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "DLPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A    -    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "DLPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A    -    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "DLPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A    -    -    -    -  C1B  C2B  C3B    -    -    -   "),
    
    "DMPC": (moltype, "    -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   -    "),
    "DMPE": (moltype, "    -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   -    "),
    
    "DNPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DNPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DNPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DNPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DNPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    
    "DOPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  D2B  C3B  C4B    -    -   "),
    "DOPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  D2B  C3B  C4B    -    -   "),
    "DOPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  D2B  C3B  C4B    -    -   "),
    "DOPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  D2B  C3B  C4B    -    -   "),
    "DOPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  D2B  C3B  C4B    -    -   "),
    
    "DPPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DPPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DPPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DPPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "DPPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "DRPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  D1B  D2B  D3B  D4B  D5B  D6B   "),
    "DRPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  D1B  D2B  D3B  D4B  D5B  D6B   "),
    "DRPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  D1B  D2B  D3B  D4B  D5B  D6B   "),
    "DRPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  D1B  D2B  D3B  D4B  D5B  D6B   "),
    "DRPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  D1B  D2B  D3B  D4B  D5B  D6B   "),

    "DSPC": (moltype, "    -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  -    "),
    "DSPE": (moltype, "    -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  -    "),

    "DTPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A    -    -    -    -    -  C1B  C2B    -    -    -    -   "),
    "DTPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A    -    -    -    -    -  C1B  C2B    -    -    -    -   "),
    "DTPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A    -    -    -    -    -  C1B  C2B    -    -    -    -   "),
    "DTPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A    -    -    -    -    -  C1B  C2B    -    -    -    -   "),
    "DTPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A    -    -    -    -    -  C1B  C2B    -    -    -    -   "),
    
    "DUPC": (moltype, "    -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   -    "),
    "DUPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  D1B  D2B  D3B  D4B  D5B    -   "),
    "DUPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  D1B  D2B  D3B  D4B  D5B    -   "),
    
    "DVPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DVPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DVPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DVPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    "DVPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  D3B  C4B    -    -   "),
    
    "DXPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DXPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DXPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DXPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    "DXPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -  C6A    -  C1B  C2B  C3B  C4B    -  C6B   "),
    
    "DYPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  C3A    -    -    -    -  C1B  D2B  C3B    -    -    -   "),
    "DYPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  C3A    -    -    -    -  C1B  D2B  C3B    -    -    -   "),
    "DYPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  C3A    -    -    -    -  C1B  D2B  C3B    -    -    -   "),
    "DYPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  C3A    -    -    -    -  C1B  D2B  C3B    -    -    -   "),
    "DYPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  C3A    -    -    -    -  C1B  D2B  C3B    -    -    -   "),

    "LPPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "LPPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "LPPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "LPPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B    -    -    -   "),
    "LPPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  C3A  C4A    -    -    -  C1B  C2B  C3B    -    -    -   "),

    "PAPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PAPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PAPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PAPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PAPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "PEPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PGPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PGPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PGPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PGPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PGPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  C2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "PIPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PIPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PIPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PIPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PIPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "POPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "POPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "POPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "POPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "POPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  C3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),

    "PQPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PQPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  C1A  D2A  D3A  C4A    -    -    -  C1B  C2B  C3B  C4B    -    -   "),
    
    "PRPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  C1B  C2B  C3B  C4B    -    -   "),
    "PRPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  C1B  C2B  C3B  C4B    -    -   "),
    "PRPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  C1B  C2B  C3B  C4B    -    -   "),
    "PRPG": (moltype, "      -    -  GL0    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  C1B  C2B  C3B  C4B    -    -   "),
    "PRPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A  D6A    -  C1B  C2B  C3B  C4B    -    -   "),

    "PUPA": (moltype, "      -    -    -    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PUPC": (moltype, "      -    -  NC3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PUPE": (moltype, "      -    -  NH3    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  C1B  C2B  C3B  C4B    -    -   "),
    "PUPS": (moltype, "      -    -  CNO    -  PO4  GL1  GL2  D1A  D2A  D3A  D4A  D5A    -    -  C1B  C2B  C3B  C4B    -    -   "),

})

# PI45P2
moltype = "PI45P2"
lipidsx[moltype] = (    0, .5,  0,  0, .5,  0,  0, .5,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0)
lipidsy[moltype] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (    7,  8,  8,  6,  9,  9,  5,  5,  4,  3,  2,  1,  0,  4,  3,  2,  1,  0)
lipidsa.update({      # 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18

})

#Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--26-26-27-28-29
#  |/   |/  |/  |/    |
#  11   8   5   2    19--20-21-22-23-24 
moltype = "GLYCOLIPIDS"
lipidsx[moltype] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1)
lipidsy[moltype] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (    6,   7,   7,   8,   9,  9, 10, 11, 11,   12,   13,   13,   10,    9,  10,    8,   11,    5,    5,   4,   3,   2,   1,   0,   4,   3,   2,   1,   0)
lipidsa.update({      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29
    "GM1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17  GM18  GM19 GM20 GM21 GM22 GM23 GM24 GM25 GM26 GM27 GM28   - "), 

    "CER" : (moltype, "  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),

    "PI"  : (moltype, " C1   C2   C3    -   CP   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   CU1  CU2  CU3  CU4  CU5"),

#lipids for thylakoid membrane of cyanobacteria: oleoyl tail at sn1 and palmiotyl chain at sn2. SQDG no double bonds

#lipids for thylakoid membrane of spinach: for the *T both chains are triple unsaturated and the *G have a triple unsaturated chain at sn1 and a palmitoyl chain at sn2. 

})

# Prototopology for mycolic acid(s)
#
#  1--2--3--4--5--6--7--8
#                       |
# 16-15-14-13-12-11-10--9
# |
# 17-18-19-20-21-22-23-24
#                     /
# 32-31-30-29-28-27-25-26
#


# Sterols
moltype = "sterol"
lipidsx[moltype] = (      0,   1,   0,   0,  1,   0, 0.5, 0.5,   0,  0)
lipidsy[moltype] = (      0,   0,   0,   0,  0,   0,   0,   0,   0,  0)
lipidsz[moltype] = (    5.3, 4.5, 3.9, 3.3,  3, 2.6, 4.5, 2.6, 1.4,  0)
lipidsa.update({
    "CHOL": (moltype, " ROH  R1  R2  R3  R4   -  R5  R6  C1  C2 "),
})

# Lists for automatic charge determination
charges = {"ARG":1, "LYS":1, "ASP":-1, "GLU":-1, "DOPG":-1, "POPG":-1, "DOPS":-1, "POPS":-1, "DSSQ":-1}

a,  b  = math.sqrt(2)/20, math.sqrt(2)/60
ct, st = math.cos(math.pi*109.47/180), math.sin(math.pi*109.47/180) # Tetrahedral

# Get a set of coordinates for a solvent particle with a given name
# Dictionary of solvents; First only those with multiple atoms
solventParticles = {
    "PW":       (("W",(-0.07,0,0)),                          # Polarizable water
                 ("WP",(0.07,0,0)),
                 ("WM",(0.07,0,0))),
    "IM":       (("SI1",(0,0,0)),                             # Imidazole
                 ("SI2",(0.2,0,0)),
                 ("SI3",(0.1,0.1,0))),
    "BMW":      (("C",(0,0,0)),
                 ("Q1",(0.12,0,0)),
                 ("Q2",(-0.06,math.cos(math.pi/6)*0.12,0))), # BMW water
    "SPC":      (("OW",(0,0,0)),                             # SPC
                 ("HW1",(0.01,0,0)),
                 ("HW2",(0.01*ct,0.01*st,0))),
    "SPCM":     (("OW",(0,0,0)),                             # Multiscale/Martini SPC 
                 ("HW1",(0.01,0,0)),
                 ("HW2",(0.01*ct,0.01*st,0)),
                 ("vW",(0,0,0))),
    "FG4W":     (("OW1",(a,a,a)),                            # Bundled water
                 ("HW11",(a,a-b,a-b)),
                 ("HW12",(a,a+b,a+b)),
                 ("OW2",(a,-a,-a)),
                 ("HW21",(a-b,-a,-a+b)),
                 ("HW22",(a+b,-a,-a-b)),
                 ("OW3",(-a,-a,a)),
                 ("HW31",(-a,-a+b,a-b)),
                 ("HW32",(-a,-a-b,a+b)),
                 ("OW4",(-a,a,-a)),
                 ("HW41",(-a+b,a,-a+b)),
                 ("HW42",(-a-b,a,-a-b))),
    "FG4W-MS":  (("OW1",(a,a,a)),                            # Bundled water, multiscaled
                 ("HW11",(a,a-b,a-b)),
                 ("HW12",(a,a+b,a+b)),
                 ("OW2",(a,-a,-a)),
                 ("HW21",(a-b,-a,-a+b)),
                 ("HW22",(a+b,-a,-a-b)),
                 ("OW3",(-a,-a,a)),
                 ("HW31",(-a,-a+b,a-b)),
                 ("HW32",(-a,-a-b,a+b)),
                 ("OW4",(-a,a,-a)),
                 ("HW41",(-a+b,a,-a+b)),
                 ("HW42",(-a-b,a,-a-b)),
                 ("VZ",(0,0,0))),
    "GLUC":     (("B1",(-0.11, 0,   0)),
                 ("B2",( 0.05, 0.16,0)),
                 ("B3",( 0.05,-0.16,0))),
    "FRUC":     (("B1",(-0.11, 0,   0)),
                 ("B2",( 0.05, 0.16,0)),
                 ("B3",( 0.05,-0.16,0))),
    "SUCR":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "MALT":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "CELL":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "KOJI":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "SOPH":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "NIGE":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "LAMI":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "TREH":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    }

# Update the solvents dictionary with single atom ones
for s in ["W","NA","CL","Mg","K","BUT","G"]:
    solventParticles[s] = ((s,(0,0,0)),)

# Apolar amino acids nd stuff for orienting proteins in membrane 
apolar = "ALA CYS PHE ILE LEU MET VAL TRP PLM CLR".split()

# molar mass of the gas
molmass={"CO2":44,"N2":28,"O2":32,"H2":2,"AIR":29}

## PRIVATE PARTS FROM THIS POINT ON ##

S = str
F = float
I = int
R = random.random

# create three-dimensional coordinates with uniform distribution on the sphere 
def balance(r,a):
    q = 34
    s = 1/(math.sqrt(5)*q)
    pi = (math.sqrt(5)-1)/2
    coord = []
    n= int(4*math.pi*pow(r,2)/a)

    for i in range(1, n):
        z = r*(2*i-1)/(n-1)
        az = abs(z-r)
        y = math.sqrt(r**2 - az**2)*math.cos(2*math.pi*i*pi)
        x = math.sqrt(r**2 - az**2)*math.sin(2*math.pi*i*pi)
        z = z -r
        coord.append([x, y, z])
    return coord

# Convert cartesian coordinates to spherical coordinates
def trans1(a, r): 
    coordi = []
    for i in a:
        x = i[0]
        y = i[1]
        z = i[2]
        theta = math.acos(z/r)
        phi = math.atan(y/x)
        if x > 0 and y > 0:
            phi = phi
        if x < 0 and y > 0:
            phi = phi + math.pi
        if x < 0 and y < 0:
            phi = phi + math.pi
        if x > 0 and y < 0:
            phi = phi + 2*math.pi
        coordi.append([phi, theta])
    return coordi

# Spherical pole projection: turn spherical coordinates to plane cartesian coordinates
def trans2(r, a):
    i = (-1)**0.5
    c = []
    for j in a:
        phi = j[0]
        theta = j[1]
        z = r*1/math.tan(theta/2)*cmath.exp(i*phi)
        x = z.real
        y = z.imag
        c.append([x, y])
    return c

# Reverse projection of plane cartesian coordinates into spherical coordinates
def trans3(z, r):
    zz = cmath.polar(z)
    zr = zz[0]
    phi = zz[1]
    theta = 2*math.atan(r/zr)
    return(phi,theta)

# Create a plural
def trans4(x,y):
    z = complex(x,y)
    return z

# Spherical coordinates to cartesian coordinates
def trans5(a,r):

    phi = a[0]
    theta = a[1]

    x = r*math.sin(theta)*math.cos(phi)
    y = r*math.sin(theta)*math.sin(phi)
    z = r*math.cos(theta)

    c = (x,y,z)

    return c

# the shortest distance between the two points
def distance(a,b):
    d = math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
    return d

def vector(v):
    if type(v) == str and "," in v:
        return [float(i) for i in v.split(",")]
    return float(v)

def vvadd(a,b):    
    if type(b) in (int,float):
        return [i+b for i in a]
    return [i+j for i,j in zip(a,b)]

def vvsub(a,b):
    if type(b) in (int,float):
        return [i-b for i in a]
    return [i-j for i,j in zip(a,b)]

def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")

def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>   atom name,   res name,     res id, chain,       x,            y,             z       
    return (S(a[12:16]),S(a[17:20]),I(a[22:26]),a[21],F(a[30:38])/10,F(a[38:46])/10,F(a[46:54])/10)

d2r = 3.14159265358979323846264338327950288/180
def pdbBoxRead(a):
    # Convert a PDB CRYST1 entry to a lattice definition.
    # Convert from Angstrom to nanometer
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac) , math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]

def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z       
    return (S(a[10:15]), S(a[5:10]),   I(a[:5]), " ", F(a[20:28]),F(a[28:36]),F(a[36:44]))

def groBoxRead(a):    
    b = [F(i) for i in a.split()] + 6*[0] # Padding for rectangular boxes
    return b[0],b[3],b[4],b[5],b[1],b[6],b[7],b[8],b[2]

def readBox(a):
    x = [ float(i) for i in a.split(",") ] + 6*[0]
    if len(x) == 12: # PDB format
        return pdbBoxRead("CRYST1 "+" ".join([str(i) for i in x]))
    else:            # GRO format
        return x[0],x[3],x[4],x[5],x[1],x[6],x[7],x[8],x[2]

class Structure:
    def __init__(self,filename=None):
        self.title   = ""
        self.atoms   = []
        self.coord   = []
        self.rest    = []
        self.box     = []        
        self._center = None

        if filename:
            lines = open(filename).readlines()
            # Try extracting PDB atom/hetatm definitions
            self.rest   = []
            self.atoms  = [pdbAtom(i) for i in lines if isPDBAtom(i) or self.rest.append(i)]
            if self.atoms:             
                # This must be a PDB file
                self.title = "THIS IS INSANE!\n"
                for i in self.rest:
                    if i.startswith("TITLE"):
                        self.title = i
                self.box   = [0,0,0,0,0,0,0,0,0]
                for i in self.rest:
                    if i.startswith("CRYST1"):
                        self.box = pdbBoxRead(i)                
            else:
                # This should be a GRO file
                self.atoms = [groAtom(i) for i in lines[2:-1]]
                self.rest  = [lines[0],lines[1],lines[-1]]
                self.box   = groBoxRead(lines[-1])
                self.title = lines[0]
            self.coord = [i[4:7] for i in self.atoms]
            self.center()

    def __nonzero__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self,s):
        for i in range(len(self)):
            self.coord[i] = vvadd(self.coord[i],s)
        return self

    def center(self,other=None):
        if not self._center:
            self._center = [ sum(i)/len(i) for i in zip(*self.coord)]
        if other:
            s = vvsub(other,self._center)
            for i in range(len(self)):
                self.coord[i] = vvadd(self.coord[i],s)
            self._center = other
            return s # return the shift
        return self._center

    def diam(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j+k*k for i,j,k in self.coord]))

    def diamxy(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j for i,j,k in self.coord]))

    def fun(self,fn):
        return [fn(i) for i in zip(*self.coord)]

# Mean of deviations from initial value
def meand(v):
    return sum([i-v[0] for i in v])/len(v)

# Sum of squares/crossproducts of deviations
def ssd(u,v):
    return sum([(i-u[0])*(j-v[0]) for i,j in zip(u,v)])/(len(u)-1)

# Parse a string for a lipid as given on the command line (LIPID[:NUMBER]) 
def parse_mol(x):
    l = x.split(":")
    return l[0], len(l) == 1 and 1 or float(l[1])

## MIJN EIGEN ROUTINE ##

# Quite short piece of code for diagonalizing symmetric 3x3 matrices :)

# Analytic solution for third order polynomial
def solve_p3( a, b, c ):
    Q,R,a3 = (3*b-a**2)/9.0, (-27*c+a*(9*b-2*a**2))/54.0, a/3.0
    if Q**3 + R**2:
        t,R13 = math.acos(R/math.sqrt(-Q**3))/3, 2*math.sqrt(-Q)
        u,v,w = math.cos(t), math.sin(t+math.pi/6), math.cos(t+math.pi/3)
        return R13*u-a3, -R13*v-a3, -R13*w-a3
    else:
        R13   = math.sqrt3(R)
        return 2*R13-a3, -R13-a3, -R13-a3

# Normalization of 3-vector
def normalize(a):
    f = 1.0/math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    return f*a[0],f*a[1],f*a[2]

# Eigenvectors for a symmetric 3x3 matrix:
# For symmetric matrix A the eigenvector v with root r satisfies
#   v.Aw = Av.w = rv.w = v.rw
#   v.(A-rI)w = v.Aw - v.rw = 0 for all w
# This means that for any two vectors p,q the eigenvector v follows from:
#   (A-rI)p x (A-rI)q
# The input is var(x),var(y),var(z),cov(x,y),cov(x,z),cov(y,z)
# The routine has been checked and yields proper eigenvalues/-vectors
def mijn_eigen_sym_3x3(a,d,f,b,c,e):
    a,d,f,b,c,e=1,d/a,f/a,b/a,c/a,e/a
    b2, c2, e2, df = b*b, c*c, e*e, d*f
    roots = list(solve_p3(-a-d-f, df-b2-c2-e2+a*(f+d), a*e2+d*c2+f*b2-a*df-2*b*c*e))
    roots.sort(reverse=True)
    ux, uy, uz = b*e-c*d, b*c-a*e, a*d-b*b
    u = (ux+roots[0]*c,uy+roots[0]*e,uz+roots[0]*(roots[0]-a-d))
    v = (ux+roots[1]*c,uy+roots[1]*e,uz+roots[1]*(roots[1]-a-d))
    w = u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0] # Cross product
    return normalize(u),normalize(v),normalize(w),roots

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

tm   = []
lipL = []
lipU = []
solv = []
gas  =[]

# HII edit - lipid definition, for extra lipid definitaions
usrmols  = []
usrheads = []
usrlinks = []
usrtails = []
usrLipHeadMapp = { # Define supported lipid head beads. One letter name mapped to atom name
    "C":  ('NC3'), # NC3 = Choline
    "E":  ('NH3'), # NH3 = Ethanolamine 
    "G":  ('GL0'), # GL0 = Glycerol
    "S":  ('CNO'), # CNO = Serine
    "P":  ('PO4')  # PO4 = Phosphate
    }
usrIndexToLetter = "A B C D E F G H I J K L M N".split() # For naming lipid tail beads 

# Description
desc = "The following are the functions that the program can perform:"

# Option list
options = [
#   option           type number default description
# HII edit - lipid definition (last options are for additional lipid specification)
(  """
Input/output related options
"""),
    ("-f",      Option(tm.append,   1,        None, "Input GRO or PDB file 1: Protein")),
    ("-o",      Option(str,         1,        None, "Output GRO file: Membrane with Protein")),
    ("-p",      Option(str,         1,        None, "Optional rudimentary topology file")),
    """
Periodic boundary conditions 
If -d is given, set up PBC according to -pbc such that no periodic
images are closer than the value given.  This will make the numbers
provided for lipids be interpreted as relative numbers. If -d is
omitted, those numbers are interpreted as absolute numbers, and the
PBC are set to fit the given number of lipids in.
""",
    ("-pbc",    Option(str,         1, "hexagonal", "PBC type: hexagonal, rectangular, square, cubic, optimal or keep")),
    ("-d",      Option(float,       1,        None, "Distance between periodic images (nm)")),
    ("-dz",     Option(float,       1,           0, "Z distance between periodic images (nm)")),
    ("-x",      Option(vector,      1,           0, "X dimension or first lattice vector of system (nm)")),
    ("-y",      Option(vector,      1,           0, "Y dimension or first lattice vector of system (nm)")),
    ("-z",      Option(vector,      1,           0, "Z dimension or first lattice vector of system (nm)")),
#    ("-box",    Option(readBox,     1,        None, "Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated")),
#    ("-n",      Option(str,         1,        None, "Index file --- TO BE IMPLEMENTED")),
    """
Membrane/lipid related options.  
The options -l and -u can be given multiple times. Option -u can be
used to set the lipid type and abundance for the upper leaflet. Option
-l sets the type and abundance for the lower leaflet if option -u is
also given, or for both leaflets if option -u is not given. The
meaning of the number depends on whether option -d is used to set up
PBC
""",
#    ("-l",      Option(lipL.append, 1,   None, "Lipid type and relative abundance (NAME[:#])")),
    ("-u",      Option(lipU.append, 1,   None, "Lipid type and relative abundance (NAME[:#])")),
    ("-a",      Option(float,       1,        0.60, "Area per lipid (nm*nm)")),
    ("-asym",   Option(int,         1,        None, "Membrane asymmetry (number of lipids)")),
    ("-hole",   Option(float,       1,        None, "Make a hole in the membrane with specified radius")),
    ("-rand",   Option(float,       1,         0.1, "Random kick size (maximum atom displacement)")),
    ("-r",      Option(float,       1,        0, "radius of the Membrane sphere")),
    ("-bd",     Option(float,       1,         0.3, "Bead distance unit for scaling z-coordinates (nm)")),
    """
Solvent related options.
""",
    ("-sol",    Option(solv.append, 1,        None, "Solvent type and relative abundance (NAME[:#])")),
    ("-gas",    Option(gas.append, 1,       None, "Solvent inside type and relative abundance (NAME[:#])")),
    ("-sold",   Option(float,       1,         0.5, "Solvent diameter")),
    ("-gasd",   Option(float,       1,         None, "Gas diameter")),
    ("-solr",   Option(float,       1,         0.1, "Solvent random kick")),
    ("-gasr",   Option(float,       1,         0.1, "Gas random kick")),
    ("-gden",   Option(float,       1,         None, "Gas density")),
    ("-excl",   Option(float,       1,         1.3, "Exclusion range (nm) for solvent addition relative to membrane surface")),
    """
Salt related options.
""",
    ("-salt",   Option(str,         1,        None, "Salt concentration")),
    ("-charge", Option(str,         1,      "auto", "Charge of system. Set to auto to infer from residue names")),
    """
Define additional lipid types (same format as in lipid-martini-itp-v01.py)
""",
    ("-alname",  Option(usrmols.append,         1,        None, "Additional lipid name, x4 letter")),
    ("-alhead",  Option(usrheads.append,        1,        None, "Additional lipid head specification string")),
    ("-allink",  Option(usrlinks.append,        1,        None, "Additional lipid linker specification string")),
    ("-altail",  Option(usrtails.append,        1,        None, "Additional lipid tail specification string")),
    ]
    
args = sys.argv[1:]

if '-h' in args or '--help' in args:
    print ("\n",__file__)
    print (desc) or ("\nSomeone ought to write a description for this script...\n")
    for thing in options:
        if type(thing)!=str:
            print("%10s  %s"%(thing[0],thing[1].description))
        else:
            print(thing)
    sys.exit()

# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])

# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

# Read in the structures (if any)    
tm    = [ Structure(i) for i in tm ]

absoluteNumbers = not options["-d"]

# HII edit - lipid definition
# Add specified lipid definition to insane lipid library
for name, head, link, tail in zip(usrmols,usrheads,usrlinks,usrtails):
    print( "Lipid name %s : %s - %s - %s" % (name, head, link, tail))

    moltype = "usr_"+name
    lipidsx[moltype] = []
    lipidsy[moltype] = []
    lipidsz[moltype] = []
    headArray = (head).split()
    linkArray = (link).split()
    tailsArray = (tail).split()
    lipidDefString = ""  

    if len(tailsArray) != len(linkArray):
        print ("Error, Number of tails has to equal number of linkers")
        sys.exit()

    # Find longest tail 
    maxTail = 0
    for cTail in tailsArray:
       if len(cTail) > maxTail:
           maxTail = len(cTail)
    cBeadZ = maxTail + len(headArray) # longest tail + linker (always x1) + lengths of all heads - 1 (as it starts on 0)

    # Add head beads
    for cHead in headArray:
        lipidsx[moltype].append(0)
        lipidsy[moltype].append(0)
        lipidsz[moltype].append(cBeadZ)
        cBeadZ -= 1
        lipidDefString += usrLipHeadMapp[cHead] + " "

    # Add linkers
    for i,cLinker in enumerate(linkArray):
        lipidsx[moltype].append(max(i-0.5,0))
        lipidsy[moltype].append(0)
        lipidsz[moltype].append(cBeadZ)
        lipidDefString += "GL" + str(i+1) + " "

    # Add tails 
    for i,cTail in enumerate(tailsArray):
        cBeadZ = maxTail - 1
        
        for j,cTailBead in enumerate(cTail):
            lipidsx[moltype].append(i)
            lipidsy[moltype].append(0)
            lipidsz[moltype].append(cBeadZ)
            cBeadZ -= 1
            lipidDefString += cTailBead + str(j+1) + usrIndexToLetter[i] + " "
   
    lipidsa[name] = (moltype,lipidDefString)
# End user lipid definition

# HII edit - lipid definition, had to move this one below the user lipid definitions to scale them to.
# First all X/Y coordinates of templates are centered and scaled (magic numbers!)
for i in lipidsx.keys():
    cx = (min(lipidsx[i])+max(lipidsx[i]))/2
    lipidsx[i] = [0.25*(j-cx) for j in lipidsx[i]]
    cy = (min(lipidsy[i])+max(lipidsy[i]))/2
    lipidsy[i] = [0.25*(j-cy) for j in lipidsy[i]]

# Periodic boundary conditions

# option -pbc keep really overrides everything
if options["-pbc"].value == "keep" and tm:
    options["-x"].value = tm[0].box[:3]
    options["-y"].value = tm[0].box[3:6]
    options["-z"].value = tm[0].box[6:]    

# options -x, -y, -z take precedence over automatic determination
pbcSetX = 0
if type(options["-x"].value) in (list,tuple):
    pbcSetX = options["-x"].value
elif options["-x"].value:
    pbcSetX = [options["-x"].value,0,0]

pbcSetY = 0
if type(options["-y"].value) in (list,tuple):
    pbcSetY = options["-y"].value
elif options["-y"].value:
    pbcSetY = [0,options["-y"].value,0]

pbcSetZ = 0
if type(options["-z"].value) in (list,tuple):
    pbcSetZ = options["-z"].value
elif options["-z"].value:
    pbcSetZ = [0,0,options["-z"].value]

a  = options["-a"].value

################
## I. PROTEIN ##
################

protein  = Structure()
wholecenter = []
prot     = []
shift    = [0] # Shift in x direction per protein

## A. NO PROTEIN ---
if not tm:

    resi = 0
    # Set the box -- If there is a hole, add its radius to the distance
    pbcx = pbcy = pbcz = options["-d"].value + (options["-hole"] and options["-hole"].value or 0)
    if "hexagonal".startswith(options["-pbc"].value):
        # Hexagonal prism -- y derived from x directly
        pbcy = math.sqrt(3)*pbcx/2
        pbcz = options["-dz"].value or options["-z"].value or options["-d"].value
    elif "optimal".startswith(options["-pbc"].value): 
        # Rhombic dodecahedron with hexagonal XY plane
        pbcy = math.sqrt(3)*pbcx/2
        pbcz = math.sqrt(6)*options["-d"].value/3
    if "rectangular".startswith(options["-pbc"].value): 
        pbcz = options["-dz"].value or options["-z"].value or options["-d"].value

    # Possibly override
    pbcx = pbcSetX and pbcSetX[0] or pbcx
    pbcy = pbcSetY and pbcSetY[1] or pbcy
    pbcz = pbcSetZ and pbcSetZ[2] or pbcz

## B. PROTEIN ---
else:

    for prot in tm:

        ## a. NO MEMBRANE --
        if not lipL:

            # A protein, but don't add lipids... Just solvate the protein
            # Maybe align along principal axes and then build a cell according to PBC
                        
            # Set PBC starting from diameter and adding distance
            if "cubic".startswith(options["-pbc"].value):
                pbcx = pbcy = pbcz = prot.diam()+options["-d"].value
            elif "rectangular".startswith(options["-pbc"].value):                
                pbcx, pbcy, pbcz = vvadd(vvsub(prot.fun(max),prot.fun(min)),options["-d"].value)
            else:
                # Rhombic dodecahedron
                pbcx = pbcy = prot.diam()+options["-d"].value
                pbcz = math.sqrt(2)*pbcx/2

            # Possibly override
            pbcx = pbcSetX and pbcSetX[0] or pbcx
            pbcy = pbcSetY and pbcSetY[1] or pbcy
            pbcz = pbcSetZ and pbcSetZ[2] or pbcz

            # Center coordinates in rectangular brick -- Add solvent next
            if len(tm) == 1:
                prot.center((0.5*pbcx, 0.5*pbcy, 0.5*pbcz))

            # Do not set an exclusion range for solvent
            options["-excl"].value = -1

        ## b. PROTEIN AND MEMBRANE --
        else:
        
            # Have to build a membrane around the protein. 
            # So first put the protein in properly.

            # Center the protein and store the shift
            shift = prot.center((0,0,0))

            ## 1. Orient with respect to membrane
            # Orient the protein according to the TM region, if requested
            # This doesn't actually work very well...
            if options["-orient"]:

                # Grid spacing (nm)
                d  = options["-od"].value
                pw = options["-op"].value

                # Determine grid size
                mx,my,mz = prot.fun(min)
                rx,ry,rz = prot.fun(lambda x: max(x)-min(x)+1e-8)

                # Number of grid cells
                nx,ny,nz = int(rx/d+0.5),int(ry/d+0.5),int(rz/d+0.5)

                # Initialize grids
                atom     = [[[0 for i in range(nz+2)] for j in range(ny+2)] for k in range(nx+2)]
                phobic   = [[[0 for i in range(nz+2)] for j in range(ny+2)] for k in range(nx+2)]
                surface  = []
                for i, (ix, iy, iz) in zip(prot.atoms,prot.coord):
                    if i[1] != "DUM":
                        jx,jy,jz = int(nx*(ix-mx)/rx), int(ny*(iy-my)/ry), int(nz*(iz-mz)/rz)
                        atom[jx][jy][jz]   += 1
                        phobic[jx][jy][jz] += (i[1].strip() in apolar)

                # Determine average density
                occupd = sum([bool(k) for i in atom for j in i for k in j])
                avdens = float(sum([sum(j) for i in atom for j in i]))/occupd

                #cgofile  = open('density.cgo',"w")
                #cgofile.write('[\n')
                for i in range(nx):
                    for j in range(ny):
                        for k in range(nz):
                            if atom[i][j][k] > 0.1*avdens:
                                # Check the neighbouring cells; If one of them is not occupied, count cell as surface
                                if not (atom[i-1][j][k] and atom[i+1][j][k] and
                                        atom[i][j-1][k] and atom[i][j+1][k] and
                                        atom[i][j][k-1] and atom[i][j][k+1]):
                                    sx,sy,sz = mx+rx*(i+0.5)/nx, my+ry*(j+0.5)/ny, mz+rz*(k+0.5)/nz
                                    sw       = (2.0*phobic[i][j][k]/atom[i][j][k])**pw
                                    surface.append((sx,sy,sz,sw))
                                    #cgofile.write("    7.0, %f, %f, %f, %f,\n"%(10*sx,10*sy,10*sz,0.25*sw))
                #cgofile.write(']\n')
                #cgofile.close()

                sx, sy, sz, w = zip(*surface)
                W             = 1.0/sum(w)

                # Weighted center of apolar region; has to go to (0,0,0) 
                sxm,sym,szm   = [sum(p)*W for p in zip(*[(m*i,m*j,m*k) for m,i,j,k in zip(w,sx,sy,sz)])]

                # Place apolar center at origin
                prot.center((-sxm,-sym,-szm))
                sx, sy, sz    = zip(*[(i-sxm,j-sym,k-szm) for i,j,k in zip(sx,sy,sz)])

                # Determine weighted deviations from centers 
                dx,dy,dz      = zip(*[(m*i,m*j,m*k) for m,i,j,k in zip(w,sx,sy,sz)]) 

                # Covariance matrix for surface
                xx,yy,zz,xy,yz,zx = [sum(p)*W for p in zip(*[(i*i,j*j,k*k,i*j,j*k,k*i) for i,j,k in zip(dx,dy,dz)])]
                
                # PCA: u,v,w are a rotation matrix
                (ux,uy,uz),(vx,vy,vz),(wx,wy,wz),r = mijn_eigen_sym_3x3(xx,yy,zz,xy,zx,yz)

                # Rotate the coordinates
                prot.coord = [(ux*i+uy*j+uz*k,vx*i+vy*j+vz*k,wx*i+wy*j+wz*k) for i,j,k in prot.coord]

            ## 4. Orient the protein in the xy-plane
            ## i. According to principal axes and unit cell
            if options["-rotate"].value == "princ":

                x, y, z = zip(*prot.coord)

                # The rotation matrix in the plane equals the transpose
                # of the matrix of eigenvectors from the 2x2 covariance
                # matrix of the positions.
                # For numerical stability we do
                # d_i     = x_i - x_0
                # mean(x) = x_0 + sum(d_i)/N =
                # var(x)  = sum((d_i - mean(d))**2)/(N-1)
                xy        = ssd(x,y)
                if xy != 0:
                    xx     = ssd(x,x)
                    yy     = ssd(y,y)
                    
                    # The eigenvalues are the roots of the 2nd order
                    # characteristic polynomial, with the coefficients
                    # equal to the trace and the determinant of the 
                    # matrix.
                    t,  d  = xx+yy, xx*yy - xy*xy
                    # The two eigenvectors form a 2D rotation matrix
                    # R = ((cos,sin),(-sin,cos)), which means that
                    # the second eigenvector follows directly from
                    # the first. We thus only need to determine one.
                    l1     = t/2 + math.sqrt(0.25*t*t-d)
                
                    ux, uy = l1-yy, xy
                    lu     = math.sqrt(ux*ux+uy*uy)
                    
                    ux    /=  lu
                    uy    /=  lu
                    
                    # Finally we rotate the system in the plane by 
                    # matrix multiplication with the transpose of 
                    # the matrix of eigenvectors
                    prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in zip(x,y,z)]

            ## ii. Randomly
            elif options["-rotate"].value == "random":
                ux   = math.cos(R()*2*math.pi)
                uy   = math.sqrt(1-ux*ux)
                prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in prot.coord]

            ## iii. Specifically
            elif options["-rotate"]:
                ux   = math.cos(float(options["-rotate"].value)*math.pi/180.)
                uy   = math.sin(float(options["-rotate"].value)*math.pi/180.)
                prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in prot.coord]

            ## 5. Determine the minimum and maximum x and y of the protein 
            pmin, pmax = prot.fun(min), prot.fun(max)
            prng       = (pmax[0]-pmin[0],pmax[1]-pmin[1],pmax[2]-pmin[2])
            center     = (0.5*(pmin[0]+pmax[0]),0.5*(pmin[1]+pmax[1]))

            wholecenter.append((center,prng))

            # Set the z-dimension
            pbcz  = pbcSetZ and pbcSetZ[2]
            # If it is not set, set pbcz to the dimension of the protein
            pbcz  = pbcz or prng[2]
            pbcz += options["-dz"].value or options["-d"].value or 0

            # At this point we should shift the subsequent proteins such that they end up
            # at the specified distance, in case we have a number of them to do
            # y-shift is always -ycenter
            # x-shift is -xmin+distance+xmax(current)
            xshft, yshft = shift[-1]-pmin[0]+(options["-d"].value or 0), -center[1]
            shift.append(shift[-1]+pmax[0]+(options["-d"].value or 0))

            ## 6. Set box (brick) dimensions
            pbcx = (options["-d"].value or 0) + prng[0]
            if "square".startswith(options["-pbc"].value):
                pbcy = pbcx
            elif "rectangular".startswith(options["-pbc"].value):
                pbcy = options["-d"].value + prng[1]
            else:
                # This goes for a hexagonal cell as well as for the optimal arrangement
                # The latter is hexagonal in the membrane plane anyway...
                pbcy  = math.cos(math.pi/6)*pbcx

            ## 7. Adjust PBC for hole
            # If we need to add a hole, we have to scale the system
            # The scaling depends on the type of PBC
            if options["-hole"]:
                if ("square".startswith(options["-pbc"].value) or 
                    "rectangular".startswith(options["-pbc"].value)):
                    scale = 1+options["-hole"].value/min(pbcx,pbcy)
                else:
                    area  = options["-hole"].value**2/math.cos(math.pi/6)
                    scale = 1+area/(pbcx*pbcy)
                pbcx, pbcy = scale*pbcx, scale*pbcy

            pbcx = pbcSetX and pbcSetX[0] or pbcx
            pbcy = pbcSetY and pbcSetY[1] or pbcy

            ## 2. Shift of protein relative to the membrane center
            if options["-dm"]:
                if options["-dm"].value < 0:
                    zshift = options["-dm"].value # - max(zip(*prot.coord)[2])
                else:                        
                    zshift = options["-dm"].value # - min(zip(*prot.coord)[2])
            elif not options["-center"]:
                zshift = -shift[2]
            else:
                zshift = 0

            # Now we center the system in the rectangular 
            # brick corresponding to the unit cell
            # If -center is given, also center z in plane
            prot += (0.5*pbcx, 0.5*pbcy, zshift)

            rescoord = []

        # And we collect the atoms
        protein.atoms.extend(prot.atoms)
        protein.coord.extend(rescoord)

    # Extract the parts of the protein that are in either leaflet
    prot_up,prot_lo = [],[]
    for ix,iy,iz in protein.coord:
        if   iz > 0 and iz <  2.4:
            prot_up.append((ix,iy))
        elif iz < 0 and iz > -2.4:
            prot_lo.append((ix,iy))

    # Current residue ID is set to that of the last atom
    resi = protein.atoms[-1][2]
    
atid      = len(protein)+1
molecules = []

# The box dimensions are now (likely) set.
# If a protein was given, it is positioned in the center of the
# rectangular brick.

# Set the lattice vectors
if ("rectangular".startswith(options["-pbc"].value) or
    "square".startswith(options["-pbc"].value) or
    "cubic".startswith(options["-pbc"].value)):
    box    = [[pbcx,0,0],[0,pbcy,0],[0,0,pbcz]]
elif not lipL:
    # Rhombic dodecahedron with square XY plane
    box    = [[pbcx,0,0],[0,pbcy,0],[0.5*pbcx,0.5*pbcx,pbcz]]
elif "hexagonal".startswith(options["-pbc"].value):
    box    = [[pbcx,0,0],[math.sin(math.pi/6)*pbcx,pbcy,0],[0,0,pbcz]]
else: # optimal packing; rhombic dodecahedron with hexagonal XY plane
    box    = [[pbcx,0,0],[math.sin(math.pi/6)*pbcx,pbcy,0],[pbcx/2,pbcy/3,pbcz]]

# Override lattice vectors if they were set explicitly
box[0] = pbcSetX or box[0]
box[1] = pbcSetY or box[1]
box[2] = pbcSetZ or box[2]

grobox = (box[0][0],box[1][1],box[2][2],
          box[0][1],box[0][2],box[1][0],
          box[1][2],box[2][0],box[2][1])

pbcx, pbcy, pbcz = box[0][0], box[1][1], box[2][2]

rx, ry, rz = pbcx+1e-8, pbcy+1e-8, pbcz+1e-8

#################
## 2. MEMBRANE ##
#################

membrane = Structure()

if lipU:
    if not options["-r"].value:
        r = pbcx / (2 * math.pi)
    if 2*options["-r"].value > options["-x"].value:
        r = pbcx / (2 * math.pi)
        print("Warning: the membrane radius is too large ")
    else:
        r = options["-r"].value

    # Create uniformly distributed rectangular coordinates
    grid = balance(r,a)
    gridcopy = copy.deepcopy(grid)

    # Number of lipids
    print(";lipids number: %d"%len(grid))

    # If there is a protein, mark the corresponding cells as occupied
    if protein: 
        for i in gridcopy:
            for j in wholecenter:
                g = trans1(i, r)
                center = j[0]
                newcenter = trans3(trans4(center[0],center[1]),r)

                d = distance(r, newcenter, g)
                p = j[1]

                if d < max(p[0], p[1]):
                    grid.remove(i)

    # Set the XY coordinates
    # To randomize the lipids we add a random number which is used for sorting
    random.seed()
    upperl = []
    upper = []
    for i in grid:
        j = trans2(r, trans1([i], r))
        upperl.append((random.random(), j))
    # Sort on the random number
    upperl.sort()

    for i in upperl:
        upper.append(i[1][0])
    
    # Upper leaflet
    lipU, numU = zip(*[ parse_mol(i) for i in lipU ])
    totU       = float(sum(numU))
    num_up     = [int(len(upper)*i/totU) for i in numU]
    lip_up     = [l for i,l in zip(num_up,lipU) for j in range(i)]
    leaf_up    = ( 1,zip(lip_up,upper))

    molecules  = list(zip(lipU,num_up))

    leaflet = leaf_up[0]
    leaf_lip = leaf_up[1]

    kick       = options["-rand"].value

    coord = []
    maxz=[]

    # Build the membrane
    for i in leaf_lip:
        lipid = i[0]
        pos = i[1]
        # Increase the residue number by one
        resi += 1
        # Set the random rotation for this lipid
        rangle = random.random() * math.pi
        rcos = math.cos(rangle)
        rsin = math.sin(rangle)
        # Fetch the atom list with x,y,z coordinates
        atoms = zip(lipidsa[lipid][1].split(), lipidsx[lipidsa[lipid][0]], lipidsy[lipidsa[lipid][0]],
                    lipidsz[lipidsa[lipid][0]])
        j = []
        for i in atoms:
            if i[0] != "-":
                j.append(i)
        # Only keep atoms appropriate for the lipid
        at, ax, ay, az = zip(*j)
        # The z-coordinates are spaced at 0.3 nm,
        # starting with the first bead at 0.15 nm
        az = [leaflet * (0.5 + (i - min(az))) * 0.3 for i in az]
        maxz.append(max(az))
        xx = zip(ax, ay)

        nx = []
        ny = []

        newx =[]
        newy =[]

        cx = []
        cy = []
        cz = []

        for i, j in xx:
            nx.append(rcos * i - rsin * j + pos[0] + random.random() * kick)
            ny.append(rsin * i - rcos * j + pos[1] + random.random() * kick)

        a = pos[0]
        b = pos[1]
        c = trans3(trans4(a, b), r)
        theta = c[1]
        k = 1 / (math.sqrt(2) * pow(math.sin(theta / 2), 2))

        for i in zip(nx,ny):

            x = i[0]
            y = i[1]

            x1 = k*x+(1-k)*a
            y1 = k*y +(1-k)*b

            newx.append(x1)
            newy.append(y1)

        for ix, iy, iz in zip(newx, newy, az):

            z = trans4(ix, iy)
            phi, theta = trans3(z, r)

            ix = r * math.sin(theta) * math.cos(phi)
            nz = r * math.cos(theta)
            iy = r * math.sin(theta) * math.sin(phi)

            ix += iz * math.sin(theta) * math.cos(phi)
            nz += iz * math.cos(theta)
            iy += iz * math.sin(theta) * math.sin(phi)

            cx.append(ix)
            cy.append(iy)
            cz.append(nz)

            # Add the atoms to the list
        for i in range(len(at)):
            atom = "%5d%-5s%5s%5d" % (resi, lipid, at[i], atid)
            membrane.coord.append((cx[i], cy[i], cz[i]))
            membrane.atoms.append((at[i], lipid, resi, 0, 0, 0))
            atid += 1

    sumx=0
    sumy=0
    sumz=0

    for i in range(len(membrane.coord)):
        sumx = sumx + membrane.coord[i][0]
        sumy = sumy + membrane.coord[i][1]
        sumz = sumz + membrane.coord[i][2]

    middle = [sumx/len(membrane.coord),sumy/len(membrane.coord),sumz/len(membrane.coord)]

    d1 = math.sqrt((pbcx/2-middle[0])**2)
    d2 = math.sqrt((pbcy/2-middle[1])**2)
    d3 = math.sqrt((pbcz/2-middle[2])**2)

    membrane += (d1,d2,d3)

################
## 3. SOLVENT ##
################

# Charge of the system so far

last = None
mcharge = 0
for j in membrane.atoms:
    if not j[0].strip().startswith('v') and j[1:3] != last:
        mcharge += charges.get(j[1].strip(),0)  
    last = j[1:3]

last = None
pcharge = 0
for j in protein.atoms:
    if not j[0].strip().startswith('v') and j[1:3] != last:
        pcharge += charges.get(j[1].strip(),0)  
    last = j[1:3]

mcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in membrane.atoms])])
pcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in protein.atoms if not j[0].strip().startswith('v')])])

charge  = mcharge + pcharge
plen, mlen, slen = 0, 0, 0
plen = protein and len(protein) or 0
mlen = membrane and len(membrane) or 0

print ("; NDX Solute %d %d" % (1, protein and plen or 0))
print ( "; Charge of protein: %f" % pcharge)

print("; NDX Membrane %d %d" % (1+plen, membrane and plen+mlen or 0))
print("; Charge of membrane: %f" % mcharge)
print( "; Total charge: %f" % charge)

if solv:

    # First get names and relative numbers for each solvent
    solnames, solnums = zip(*[parse_mol(i) for i in solv])
    solnames, solnums = list(solnames), list(solnums)
    totS = float(sum(solnums))

    gasnames, gasnums = zip(*[parse_mol(i) for i in gas])
    gasnames, gasnums = list(gasnames), list(gasnums)
    totG = float(sum(gasnums))

    # Set up a grid
    d        = 1/options["-sold"].value
    dds = options["-sold"].value

    nx,ny,nz = int(1+d*pbcx),int(1+d*pbcy),int(1+d*pbcz)
    dx,dy,dz = pbcx/nx,pbcy/ny,pbcz/nz
    excl,hz  = int(nz*options["-excl"].value/pbcz), int(0.5*nz)

    zshift   = 0

    if membrane:
        memz   = [i[2] for i in membrane.coord]
        midz   = (max(memz)+min(memz))/2
        hz     = int(nz*midz/pbcz)  # Grid layer in which the membrane is located
        zshift = (hz+0.5)*nz - midz # Shift of membrane middle to center of grid layer

    # Initialize a grid of solvent, spanning the whole cell
    # Exclude all cells within specified distance from membrane center
    grid   = [[[i for i in range(nz)] for j in range(ny)] for i in range(nx)]

    if gas:
        NA = 6.023 * (10 ** 23)
        # Set up a grid
        if not options["-gasd"].value and not options["-gden"].value:
            print("Warning: You have to put in the gas density or the atomic distance")
            sys.exit()

        if options["-gasd"].value:
            fddg = options["-gasd"].value

        if options["-gden"].value:
            density = options["-gden"].value
            for x in gasnames:
                fddg = pow((molmass[x] / (density * NA * (10 ** (-24)))), (1 / 3))

        if options["-gasd"].value and options["-gden"].value:
            print("Warning: You have entered two quantities of the same meaning")
            sys.exit()

        nxg, nyg, nzg = int(1 + (1 / fddg) * pbcx), int(1 + (1 / fddg) * pbcy), int(1 + (1 / fddg) * pbcz)
        dxg, dyg, dzg = pbcx / nxg, pbcy / nyg, pbcz / nzg
        hzg = int(0.5 * nzg)

        zshiftg   = 0

        if membrane:
            memz   = [i[2] for i in membrane.coord]
            midz   = (max(memz)+min(memz))/2
            hzg     = int(nzg*midz/pbcz)  # Grid layer in which the membrane is located
            zshiftg = (hzg+0.5)*nzg - midz # Shift of membrane middle to center of grid layer
    
        # Initialize a grid of solvent, spanning the whole cell
        # Exclude all cells within specified distance from membrane center
        gridg   = [[[i for i in range(nzg) ] for j in range(nyg)] for i in range(nxg)]

    # Flag all cells occupied by protein or membrane
    for x,y,z in protein.coord+membrane.coord:
        if z >= pbcz:
            x -= box[2][0]
            y -= box[2][1]
            z -= box[2][2]
        if z < 0:
            x += box[2][0]
            y += box[2][1]
            z += box[2][2]
        if y >= pbcy: 
            x -= box[1][0]
            y -= box[1][1]
        if y < 0: 
            x += box[1][0]
            y += box[1][1]
        if x >= pbcx: 
            x -= box[0][0]
        if x < 0: 
            x += box[0][0]
        grid[int(nx*x/rx)][int(ny*y/ry)][int(nz*z/rz)] = False
        gridg[int(nxg*x/rx)][int(nyg*y/ry)][int(nzg*z/rz)]=False

    # Flag all cells inside the membrane and put gas inside
    rad   = r
    outrad = max(maxz)+rad
    # grid1=[]

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if ((i*dds-pbcx/ 2)**2+(j*dds-pbcy/ 2)**2+(k*dds-pbcz/ 2)**2)<=(outrad+0.5*excl)**2:
                    grid[i][j][k]=False
                    
    for i in range(nxg):
        for j in range(nyg):
            for k in range(nzg):
                if ((i*fddg - pbcx/ 2) ** 2 + (j*fddg  - pbcy / 2 ) ** 2 + (k*fddg - pbcz/2 ) ** 2) >= rad ** 2:
                    gridg[i][j][k]=False

    # Set the center for each solvent molecule
    kick = options["-solr"].value
    grid = [(R(), (i + 0.5 + R() * kick) * dx, (j + 0.5 + R() * kick) * dy, (k + 0.5 + R() * kick) * dz)
            for i in range(nx) for j in range(ny) for k in range(nz) if grid[i][j][k] ]
    kickg = options["-gasr"].value
    gridg = [(R(), (i + 0.5 + R() * kickg) * dxg, (j + 0.5 + R() * kickg) * dyg, (k + 0.5 + R() * kickg) * dzg)
            for i in range(nxg) for j in range(nyg) for k in range(nzg) if gridg[i][j][k] ]
    
    # Sort on the random number
    grid.sort()
    gridg.sort()

    # 'grid' contains all positions on which a solvent molecule can be placed.
    # The number of positions is taken as the basis for determining the salt concentration.
    # This is fine for simple salt solutions, but may not be optimal for complex mixtures
    # (like when mixing a 1M solution of this with a 1M solution of that

    # Set the number of ions to add
    nna, ncl = 0, 0
    if options["-salt"]:

        # If the concentration is set negative, set the charge to zero
        if options["-salt"].value.startswith("-"):
            charge = 0
            options["-salt"].value = -float(options["-salt"].value)
        else:
            options["-salt"].value = float(options["-salt"].value)

        # Determine charge to use, either determined or given on command line
        if options["-charge"].value != "0":
            charge = (options["-charge"].value != "auto") and int(options["-charge"].value) or charge
        else:
            charge = 0

        # Determine number of sodium and chloride to add
        concentration = options["-salt"].value
        nsol = ("SPC" in solnames and 1 or 4)*len(grid)
        ncl  = max(max(0,charge),int(.5+.5*(concentration*nsol/(27.7+concentration)+charge)))
        nna  = ncl - charge
                        
    # Correct number of grid cells for placement of solvent
    ngrid   = len(grid) - nna - ncl
    num_sol = [int(ngrid*i/totS) for i in solnums]
    num_gas= [int(len(gridg)*i/totG) for i in gasnums]

    # Add salt to solnames and num_sol
    if nna:
        solnames.append("NA")
        num_sol.append(nna)
        solv.append("NA")
    if ncl:
        solnames.append("CL")
        num_sol.append(ncl)
        solv.append("CL")

    # Names and grid positions for solvent molecules
    solvent    = list(zip([s for i,s in zip(num_sol,solnames) for j in range(i)],grid))
    gasin      = list(zip([s for i,s in zip(num_gas,gasnames) for j in range(i)],gridg))

    # Extend the list of molecules (for the topology)
    molecules.extend(list(zip(solnames,num_sol)))
    molecules.extend(list(zip(gasnames,num_gas)))

    # Build the solvent
    sol = []
    gas = []
    for resn,(rndm,x,y,z) in solvent:
        resi += 1
        solmol = solventParticles.get(resn)
        if solmol and len(solmol) > 1:       
            # Random rotation (quaternion)
            u,  v,  w       = random.random(), 2*math.pi*random.random(), 2*math.pi*random.random()
            s,  t           = math.sqrt(1-u), math.sqrt(u)
            qw, qx, qy, qz  = s*math.sin(v), s*math.cos(v), t*math.sin(w), t*math.cos(w)
            qq              = qw*qw-qx*qx-qy*qy-qz*qz         
            for atnm,(px,py,pz) in solmol:                
                qp = 2*(qx*px + qy*py + qz*pz)
                rx = x + qp*qx + qq*px + qw*(qy*pz-qz*py)
                ry = y + qp*qy + qq*py + qw*(qz*px-qx*pz)
                rz = z + qp*qz + qq*pz + qw*(qx*py-qy*px)
                sol.append(("%5d%-5s%5s%5d"%(resi%1e5,resn,atnm,atid%1e5),(rx,ry,rz)))
                atid += 1
        else:          
            sol.append(("%5d%-5s%5s%5d"%(resi%1e5,resn,solmol and solmol[0][0] or resn,atid%1e5),(x,y,z)))
            atid += 1

    for resn,(rndm,x,y,z) in gasin:
        resi += 1
        solmol = solventParticles.get(resn)
        if solmol and len(solmol) > 1:
            # Random rotation (quaternion)
            u,  v,  w       = random.random(), 2*math.pi*random.random(), 2*math.pi*random.random()
            s,  t           = math.sqrt(1-u), math.sqrt(u)
            qw, qx, qy, qz  = s*math.sin(v), s*math.cos(v), t*math.sin(w), t*math.cos(w)
            qq              = qw*qw-qx*qx-qy*qy-qz*qz
            for atnm,(px,py,pz) in solmol:
                qp = 2*(qx*px + qy*py + qz*pz)
                rx = x + qp*qx + qq*px + qw*(qy*pz-qz*py)
                ry = y + qp*qy + qq*py + qw*(qz*px-qx*pz)
                rz = z + qp*qz + qq*pz + qw*(qx*py-qy*px)
                sol.append(("%5d%-5s%5s%5d"%(resi%1e5,resn,atnm,atid%1e5),(rx,ry,rz)))
                atid += 1
        else:
            sol.append(("%5d%-5s%5s%5d" % (resi % 1e5, resn, solmol and solmol[0][0] or resn, atid % 1e5), (x, y, z)))
            atid += 1

#Write the output
slen = solvent and len(sol) or 0
print("; NDX Solvent %d %d" % (1+plen+mlen, solvent and plen+mlen+slen or 0))
print("; NDX System %d %d" % (1, plen+mlen+slen))

# Open the output stream
oStream = options["-o"] and open(options["-o"].value,"w") or sys.stdout

# Print the title
if membrane.atoms:
    title  = "INSANE! Membrane UpperLeaflet>"+":".join(lipU)+"="+":".join([str(i) for i in numU])
    # title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in numL])

    if protein:
        title = "Protein in " + title
else:
    title = "Insanely solvated protein."
print(title,file = oStream)

# Print the number of atoms
print("%5d"%(len(membrane)+len(sol)),file = oStream)

# Print the atoms
id = 1
if protein:
    for i in range(len(protein)):
        at,rn,ri = protein.atoms[i][:3]
        x,y,z    = protein.coord[i]
        oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id,x,y,z))
        id += 1
if membrane:
    for i in range(len(membrane)):
        at,rn,ri = membrane.atoms[i][:3]
        x,y,z    = membrane.coord[i]
        oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id,x,y,z))
        id += 1
if sol:
    # Print the solvent
    print("\n".join([i[0]+"%8.3f%8.3f%8.3f"%i[1] for i in sol]),file = oStream)

print("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"%grobox,file = oStream)

if options["-p"].value:
     # Write a rudimentary topology file
     top = open(options["-p"].value,"w")
     print('#include "martini.itp"\n',file=top)
     print('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title,file=top)
     if protein:
         print("%-10s %5d"%("Protein",1),top)
     print("\n".join("%-10s %7d"%i for i in molecules),file=top)
     top.close()
else:
     print("\n".join("%-10s %7d"%i for i in molecules),file=sys.stderr)

# Auto-generate ndx file
if options["-o"] and options["-o"].value and options["-o"].value != sys.stdout:
    try:
        import os

        if hasattr(options["-o"].value, 'name'):
            gro_path = options["-o"].value.name
        else:
            gro_path = options["-o"].value

        ndx_path = os.path.splitext(gro_path)[0] + ".ndx"

        all_resnames = [mol[0] for mol in molecules]

        lipid_keywords = set(lipidsa.keys())
        gas_keywords = {"CO2", "N2", "O2", "H2", "AIR"}
        water_keywords = {"W", "PW", "SPC", "SPCM", "FG4W", "FG4W-MS", "BMW"}
        ion_keywords = {"NA", "CL", "Mg", "K"}

        lipid_resnames = [r for r in all_resnames if r in lipid_keywords]
        gas_resnames = [r for r in all_resnames if r in gas_keywords]
        water_resnames = [r for r in all_resnames if r in water_keywords]
        ion_resnames = [r for r in all_resnames if r in ion_keywords]

        generate_ndx(
            gro_path=gro_path,
            ndx_path=ndx_path,
            lipid_resnames=lipid_resnames,
            gas_resnames=gas_resnames,
            water_resnames=water_resnames,
            ion_resnames=ion_resnames,
        )

    except ImportError:
        print("\n; Warning: ndx_generator.py not found. Skipping index file generation.")
    except Exception as e:
        print(f"\n; Error generating index file: {e}")
else:
    if options["-o"]:
        print("\n; Warning: Output option '-o' not configured for ndx generation.")