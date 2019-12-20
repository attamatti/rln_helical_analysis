#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
import math
###---------function: read the star file get the header, labels, and data -------------#######

def read_3p1star(infile):
	inDMCloop = False
	labels,data = {},[]
	dta = open(infile,'r').readlines()
	for i in dta:
		if i.split()== ['data_particles']:
			inDMCloop = True
		elif i.strip('\n') == '# version 30001':
			inDMCloop = False
		if inDMCloop == True:
			if '_rln' in i:
				labels[i.split()[0]] = int(i.split()[1].replace('#',''))-1
			elif len(i.split()) == len(labels) and len(i.split()) > 0:		
				data.append(i.split())
	return(labels,data)


#---------------------------------------------------------------------------------------------#

def graph_psi(lines,plotname):
	psis,locs,rots = [],[],[]
	for i in lines:
		psis.append(float(i[labels['_rlnAnglePsi']]))
		rots.append(float(i[labels['_rlnAngleRot']]))
		locs.append(float(i[labels['_rlnHelicalTrackLengthAngst']]))
		
	plt.scatter(locs,psis,c='blue')
	plt.scatter(locs,rots,c='red')
	plt.savefig('{0}.png'.format(plotname))
	plt.close()

## read the starfile divide into individual micrographs
errmsg = 'USAGE: rln3p0_psi_analysis.py <starfile> <minimum number of segments>'
try:
	labels,data = read_3p1star(sys.argv[1])
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
			hpls.append(float(j[labels['_rlnHelicalTrackLengthAngst']]))
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