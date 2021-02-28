import sys
import os
import glob
import re
from prody import *
from matplotlib.pylab import *
import pandas as pd

# arguments
pdbname = sys.argv[1:]

files = glob.glob(list(pdbname)[0]+'*'+'.pdb')
files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

hinges_df = pd.DataFrame([])
for f in files:
	#print('check order: ' + str(f))
	prot = parsePDB(f)
	calphas = prot.select('calpha')
	
	gnm = GNM('prot')
	gnm.buildKirchhoff(calphas)
	gnm.calcModes(1)
	
	slowest_hinges = gnm[0].getHinges()
	
	hinges = {'Resname': calphas[slowest_hinges].getResnames(),
			'Resid': calphas[slowest_hinges].getResnums()}
	hinges = pd.DataFrame(hinges, columns = ['Resname', 'Resid'])
	hinges = hinges['Resname'].str.capitalize() + hinges['Resid'].astype(str)
	
	hinges_df = hinges_df.append(hinges, ignore_index=True)
	#print(hinges_df)

hinges_df.index = files
hinges_df.to_csv(list(pdbname)[0]+'_hinges'+'.csv', index = True, header = False)