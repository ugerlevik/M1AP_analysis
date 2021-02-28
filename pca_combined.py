import sys
from prody import *
from matplotlib.pylab import *

import matplotlib.pyplot as plt

case, num  = sys.argv[1:]

str1 = parsePDB("../3_WT/3_production/seokWT_onlyProt.pdb")
ens1 = Trajectory("../3_WT/3_production/seokWT_250ns_onlyProt_each25k.dcd")
ens1.addFile("../3_WT/4_production_repeat2/seokWT_250ns_r2_onlyProt_each25k.dcd")

str2 = parsePDB("../"+num+"_"+case+"/3_production/"+"seok"+case+"_onlyProt.pdb")
ens2 = Trajectory("../"+num+"_"+case+"/3_production/"+"seok"+case+"_250ns_onlyProt_each25k.dcd")
ens2.addFile("../"+num+"_"+case+"/4_production_repeat2/"+"seok"+case+"_250ns_r2_onlyProt_each25k.dcd")

ens1.link(str1)
ens1.setCoords(str1)
ens1.setAtoms(str1.calpha)
#ens1.superpose()

ens2.link(str2)
ens2.setCoords(str2)
ens2.setAtoms(str2.calpha)
#ens2.superpose()

eda_ens1 = EDA("WT"+' Ensemble')
eda_ens1.buildCovariance( ens1 )
eda_ens1.calcModes()

eda_ens2 = EDA(case+' Ensemble')
eda_ens2.buildCovariance( ens2 )
eda_ens2.calcModes()

showProjection(ens1, eda_ens1[:2], rmsd = False, color = 'black', marker = '.', label = 'Wild-type', alpha = 0.5)
showProjection(ens2, eda_ens2[:2], rmsd = False, color = 'red', marker = '.', label = case, alpha = 0.5)

plt.legend(loc='upper left')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.xlim(-175,200)
plt.ylim(-175,300)
plt.savefig('../_outputs/PCA_combined/PCA_WT_'+case+'_2D.pdf', bbox_inches = 'tight', pad_inches = 0)
