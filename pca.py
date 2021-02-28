import sys
from prody import *
from matplotlib.pylab import *

import matplotlib.pyplot as plt

case, num  = sys.argv[1:]

str1 = parsePDB("../3_WT/3_production/seokWT_onlyProt.pdb")
ens1 = parseDCD("../3_WT/3_production/seokWT_250ns_onlyProt_each25k.dcd")

str2 = parsePDB("../"+num+"_"+case+"/3_production/"+"seok"+case+"_onlyProt.pdb")
ens2 = parseDCD("../"+num+"_"+case+"/3_production/"+"seok"+case+"_250ns_onlyProt_each25k.dcd")

ens1.setCoords(str1)
ens1.setAtoms(str1.calpha)
ens1.superpose()

ens2.setCoords(str2)
ens2.setAtoms(str2.calpha)
ens2.superpose()

eda_ens1 = EDA("WT"+' Ensemble')
eda_ens1.buildCovariance( ens1 )
eda_ens1.calcModes()

eda_ens2 = EDA(case+' Ensemble')
eda_ens2.buildCovariance( ens2 )
eda_ens2.calcModes()

# showProjection(ens1, eda_ens1[:3], color = 'black', marker = '.', label = 'Wild-type', alpha = 0.5)
# showProjection(ens2, eda_ens2[:3], color = 'red', marker = '.', label = id2, alpha = 0.5)
# plt.legend(loc='upper left')
# plt.savefig('PCA_'+id1+'_'+id2+'_3D.pdf')
# plt.close()

showProjection(ens1, eda_ens1[:2], rmsd = False, color = 'black', marker = '.', label = 'Wild-type', alpha = 0.5)
showProjection(ens2, eda_ens2[:2], rmsd = False, color = 'red', marker = '.', label = case, alpha = 0.5)
plt.legend(loc='upper right')
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.xlim(-40,40)
#plt.ylim(-30,40)
plt.savefig('../_outputs/PCA/PCA_WT_'+case+'_2D.pdf', bbox_inches = 'tight', pad_inches = 0)
