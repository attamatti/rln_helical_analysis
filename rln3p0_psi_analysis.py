#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
import math
###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile_new(f):
    inhead = True
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    count = 0
    labcount = 0
    for i in alldata:
        if '_rln' in i:
            labelsdic[i.split()[0]] = labcount
            labcount +=1
        if inhead == True:
            header.append(i.strip("\n"))
            if '_rln' in i and '#' in i and  '_rln' not in alldata[count+1] and '#' not in alldata[count+1]:
                inhead = False
        elif len(i.split())>=1:
            data.append(i.split())
        count +=1
    
    return(labelsdic,data)
#---------------------------------------------------------------------------------------------#

def graph_psi(lines,plotname):
	psis,locs,rots = [],[],[]
	for i in lines:
		psis.append(float(i[labels['_rlnAnglePsi']]))
		rots.append(float(i[labels['_rlnAngleRot']]))
		locs.append(float(i[labels['_rlnHelicalTrackLength']]))
		
	plt.scatter(locs,psis,c='blue')
	plt.scatter(locs,rots,c='red')
	plt.savefig('{0}.png'.format(plotname))
	plt.close()

## read the starfile divide into individual micrographs
errmsg = 'USAGE: rln3p0_psi_analysis.py <starfile> <minimum number of segments>'
try:
	labels,data = read_starfile_new(sys.argv[1])
except:
	sys.exit('\nERROR: Failed to read star file\n{0}'.format(errmsg))
try:
	minnosegs = int(sys.argv[2])
except:
	sys.exit('\nERROR: invalid minimum number of segemnts\n{0}'.format(errmsg))
micrographs = {}		#{micrograph:[dataline,dataline,dataline]}
for i in data:
	micrograph = i[labels['_rlnMicrographName']]
	try:
		micrographs[micrograph].append(i)
	except:
		micrographs[micrograph] = [i]
print('analyzing {0} micrographs'.format(len(micrographs)))

### separate each mcirograph into individual filaments
## iterate over all the filaments to get the maximum length for compariable graphs
hpls = []
filcount,gttencount=0,0

for micrograph in micrographs:
	mgname = micrograph.split('/')[-1].split('.')[0]
	filaments = {}				#{filamentno:[dataline,dataline,dataline]
	for i in micrographs[micrograph]:
		fil = i[labels['_rlnHelicalTubeID']]
		try:
			filaments[fil].append(i)
		except:
			filaments[fil] = [i]
			filcount+=1
	for i in filaments:
		for j in filaments[i]:
			hpls.append(float(j[labels['_rlnHelicalTrackLength']]))
maxlength = max(hpls)
print('{0} total filaments'.format(filcount))

for micrograph in micrographs:
	mgname = micrograph.split('/')[-1].split('.')[0]
	filaments = {}				#{filamentno:[dataline,dataline,dataline]
	for i in micrographs[micrograph]:
		fil = i[labels['_rlnHelicalTubeID']]
		try:
			filaments[fil].append(i)
		except:
			filaments[fil] = [i]
	for i in filaments:			
		if len(filaments[i]) > minnosegs:
			print('filament:',i,'segments:',len(filaments[i]))
			graph_psi(filaments[i],'{0}-{1}_psi'.format(mgname,i))
