# 用来储存常量
ALL_LIPIDS = ["DPPC",
     "DHPC",
     "DLPC",
     "DMPC",
     "DSPC",
     "D6PC",
     "PSM",
     "POPC",
     "POP5",
     "DOPC",
     "DAPC",
     "DDPC",
     "DUPC",
     "PDPC",
     "ODPC",
     "DEPC",
     "DPPE",
     "DHPE",
     "DLPE",
     "DMPE",
     "DSPE",
     "POPE",
     "DOPE",
     "PPCS",
     "DOPG",
     "POPG",
     "DOPS",
     "POPS",
     "DPPS",
     "DPPG",
     "CPG",
     "PPG",
     "PPT",
     "DSMG",
     "DSDG",
     "DSSQ",
     "GMO",
     ]

# Glycolipds
GLYCOLIPIDS = [
     "GM1",
     "DGDG",
     "MGDG",
     "SQDG",
     "CER",
     "GCER",
     "DPPI",
     "PI",
     "PI34",
     "CDGDG",
     "CMGDG",
     "CSQDG",
     "CSQDB",
     "PDGDG",
     "PDGDT",
     "PMGDG",
     "PMGDT",
     "PSQDG",
     ]

# Sterols
STEROLS = ['CHOL','CHOL0','CHOL1','CHOL2']

ALL_P = [*ALL_LIPIDS, *GLYCOLIPIDS, *STEROLS]
# Lipids
PC_LIPIDS = [i for i in ALL_LIPIDS if (i[-2:] == 'PC')]
PE_LIPIDS = [i for i in ALL_LIPIDS if (i[-2:] == 'PE')]
PG_LIPIDS = [i for i in ALL_LIPIDS if (i[-2:] == 'PG')]
PS_LIPIDS = [i for i in ALL_LIPIDS if (i[-2:] == 'PS')]
OTHER_LIPIDS = list(set(ALL_LIPIDS) - set(PC_LIPIDS+PE_LIPIDS+PG_LIPIDS+PS_LIPIDS))

lipids = {'PC Lipids': PC_LIPIDS,
          'PE Lipids': PE_LIPIDS,
          'PG Lipids': PG_LIPIDS,
          'PS Lipids': PS_LIPIDS,
          'Other Lipids': OTHER_LIPIDS,
          'Sterol': STEROLS,
          'Glycoipids': GLYCOLIPIDS
          }
