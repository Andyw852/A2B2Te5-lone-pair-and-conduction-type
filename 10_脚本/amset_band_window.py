"""Which bands AMSET keeps for a given energy_cutoff, and whether they line up with an existing deformation.h5.
Run: python amset_band_window.py DELIVERY [cutoff ...]      (default cutoffs 1.0 1.5)
Rule: a band is kept if any of its energies lies in [VBM - e, CBM + e] (unscissored DFT edges); it reproduces the
'Including bands 59-78' (Pb2Sb2Te5) and '61-80' (Sn2Sb2Te5) of the dense amset wave runs at e = 1.0 eV.
AMSET only checks the band COUNT of deformation.h5, so equal counts with shifted indices would pass silently; this compares indices.
"""
from pathlib import Path
import sys,json
import numpy as np
import h5py
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm

def main(D,cuts):
 for m in sm.M:
  u=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform';nel,K,E=sm.read_eigenval(u/'EIGENVAL')
  vbm,cbm=E[:,:nel].max(),E[:,nel:].min();d=D/'05_SOC_AMSET补算'/m/'03_SOC形变势'
  be=json.loads((d/'band_edges.json').read_text());h=be['hole']['hits'][0];off=h['band_global']-h['band_h5']
  nb=h5py.File(d/'deformation.h5','r')['deformation_potentials_up'].shape[0];old=(off+1,off+nb)
  for e in cuts:
   sel=np.where(((E>vbm-e)&(E<cbm+e)).any(0))[0];new=(int(sel[0])+1,int(sel[-1])+1)
   print(f'{m}  e={e:.2f} eV  dense bands {new[0]}-{new[1]} ({len(sel)})  deformation.h5 {old[0]}-{old[1]} ({nb})  '+('MATCH' if new==old else 'differ'))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),[float(x) for x in sys.argv[2:]] or [1.0,1.5])
