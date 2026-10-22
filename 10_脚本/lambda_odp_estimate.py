"""Order-of-magnitude check (19_ 1-7): is the optical deformation potential of the B-site off-centring (A1g) coordinate worth an AMSET ODP run?
Run: python lambda_odp_estimate.py DELIVERY
Deformation: band energies at fixed k between the baseline (lambda = 1) and the 15_ lambda < 1 EIGENVAL, for every state within 50 meV of
each edge, divided by the B-site displacement. Eigenvalues are referenced to the average electrostatic potential (VASP), as in frozen-phonon ODP.
Mode: only the two B atoms move (+u/-u), stiffness k from E(0) - E(1) = k u^2 (15_ section 5), omega = sqrt(k / M_B), partial density rho_B = 2 M_B / V.
Rate: nonpolar optical deformation-potential scattering of a parabolic valley at E = kT (Lundstrom), with the largest |dE/du| found and m_d = 0.5 m_e
(both upper bounds), compared with the total AMSET rate 1/tau = e / (mu m_c) at 300 K, 1e19.
"""
from pathlib import Path
import sys,csv,json
import numpy as np
import ase.io
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm

OUT='15_受控畸变/B位偏心模式_光学形变势估计.csv'
HB,E0,ME,AMU,KB=1.054571817e-34,1.602176634e-19,9.1093837e-31,1.66053907e-27,8.617333e-5
DE_RECENTRE={'Pb2Sb2Te5':0.1303,'Pb2Bi2Te5':0.1615,'Sn2Sb2Te5':0.1473,'Sn2Bi2Te5':0.1798}   # E(0) - E(1), eV/cell, 15_ section 5
MB={'Sb':121.760,'Bi':208.980}
MD=0.5

def odp_rate(D,omega,rho,T=300.):
 hw=HB*omega;E=KB*T*E0;N=1/np.expm1(hw/(KB*T*E0))
 return (D*E0*1e10)**2*(MD*ME)**1.5/(np.sqrt(2)*np.pi*HB**3*rho*omega)*(N*np.sqrt(E+hw)+(N+1)*np.sqrt(max(E-hw,0)))

def main(D):
 rows=[];mc={x['material']:x for x in csv.DictReader(open(D/'14_散射机制分项迁移率/散射输入参数.csv',encoding='utf-8-sig'))}
 off={(x['material'],x['structure'],x['cation']):float(x['offcentre_A']) for x in csv.DictReader(open(D/'15_受控畸变/键长与偏心.csv',encoding='utf-8-sig'))}
 for m,lam in (('Pb2Sb2Te5','0.50'),('Pb2Sb2Te5','0.00'),('Pb2Bi2Te5','0.00'),('Sn2Sb2Te5','0.00'),('Sn2Bi2Te5','0.00')):
  b='Bi' if 'Bi' in m else 'Sb';base=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform'
  nel,K1,E1=sm.read_eigenval(base/'EIGENVAL');_,K0,E0_=sm.read_eigenval(D/'15_受控畸变'/m/f'lam{lam}'/'EIGENVAL.gz')
  assert np.allclose(K1,K0)
  u1=off[(m,'lam1.00',b)];du=u1-off[(m,f'lam{lam}',b)]
  k=DE_RECENTRE[m]/u1**2*E0/1e-20;omega=np.sqrt(k/(MB[b]*AMU))
  rho=2*MB[b]*AMU/(ase.io.read(base/'POSCAR',format='vasp').get_volume()*1e-30)
  am=json.load(open(next((D/'14_散射机制分项迁移率/密网格AMSET'/m).glob('transport*.json'))))
  r={'material':m,'lambda':lam,'du_A':round(du,4),'omega_B_THz':round(float(omega/2/np.pi/1e12),2)}
  for edge,car,sg,band,ref in (('VB','hole',1,nel-1,lambda e:e>e.max()-0.05),('CB','electron',-1,nel,lambda e:e<e.min()+0.05)):
   i=np.where(ref(E1[:,band]))[0];s=np.abs(E0_[i,band]-E1[i,band])/du
   j=int(np.argmin(np.abs(np.array(am['doping'])-sg*1e19)));mu=np.array(am['mobility']['overall'])[j,0,0,0]*1e-4
   tot=E0/(mu*float(mc[m][f'm_c_inplane_{car}_1e19'])*ME)
   r.update({f'{edge}_mean_dEdu_eV_A':round(float(s.mean()),2),f'{edge}_max_dEdu_eV_A':round(float(s.max()),2),
             f'{car}_ODP_over_total_rate':float('%.2g'%(odp_rate(s.max(),omega,rho)/tot))})
  rows.append(r);print(r)
 with open(D/OUT,'w',newline='',encoding='utf-8-sig') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
