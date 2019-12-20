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

def serach_sawtooth_range(low,high,period,ys,vals,amplitude,interval,sign):
	difdic= {}				# {phase value: sum diffence}
	xsrange = range(0,period)
	for p in range(int(low*period),int(high*period)+10,interval):
		for n in xsrange:
			funct = '{0}*(2*amplitude/math.pi)*np.arctan((np.cos((x+n)*math.pi/p)/np.sin((x+n)*math.pi/p)))'.format(sign)
			x = np.array(vals)
			calcy = np.array(eval(funct))
			diff = sum(abs(np.array(ys) - calcy))
			difdic['{0}-{1}'.format(n,p)] = diff/len(calcy)
	return(difdic)

def draw_sawtooth(amplitude,period,vals,phase,sign):
    '''draw the sawtooth function on the graph'''
    funct = '{0}*(2*amplitude/math.pi)*np.arctan((np.cos((x+phase)*math.pi/period)/np.sin((x+phase)*math.pi/period)))'.format(sign)
    x = np.array(vals)
    y = eval(funct)
    plt.plot(x,y)
def eval_sawtooth(amplitude,period,vals,ys):
	'''evaluate the sawtooth functions and get the base phase and period values'''
	diffs= {}				# {phase value: sum diffence}
	## initial search 15% +- the expected period
	print('initial search +- 15% of {0}'.format(period))
	diffs = serach_sawtooth_range(0.85,1.15,period,ys,vals,amplitude,10,1)
	phaseperiod = min(diffs, key=diffs.get)
	####if outside range search up or down
	n=1
	inc=0.15
	while float(phaseperiod.split('-')[1]) == (1.0+(n*inc))*period:
		n+=1
		print('initial search > +{0}% range incrementing by +15%'.format(n*inc*100))
		#print('range = {},{}'.format(period*(1+((n-1)*inc)),period*(1+n*inc)))
		diffs = serach_sawtooth_range(1+((n-1)*inc),1+n*inc,period,ys,vals,amplitude,10,1)
		phaseperiod = min(diffs, key=diffs.get)

	while float(phaseperiod.split('-')[1]) == (1.0-(n*inc))*period:
		n+=1
		print('initial search < -{0}% range incrementing by -15%'.format(n*inc*100))
		diffs = serach_sawtooth_range(1-n*inc,1-((n-1)*inc),period,ys,vals,amplitude,10,1)
		phaseperiod = min(diffs, key=diffs.get)
	
	print('initial [phase,period],fitscore',phaseperiod.split('-'),diffs[phaseperiod])
	print('optimizing fit')
	diffs = serach_sawtooth_range(0.90,1.10,int(phaseperiod.split('-')[1]),ys,vals,amplitude,2,1)
	phaseperiod = min(diffs, key=diffs.get)
	
	### print top 10 hits for debugging
	dlist = list(diffs)
	dlist.sort(key=lambda x: diffs[x])
	for i in dlist[0:11]:
		print('phase: {0:03f} period: {1} score: {2}'.format(float(i.split('-')[0]),int(i.split('-')[1]),diffs[i]))

	#### do it the opposite direction

	diffs2= {}				# {phase value: sum diffence}
	## initial search 15% +- the expected period
	print('initial search +- 15% of {0}'.format(period))
	diffs2 = serach_sawtooth_range(0.85,1.15,period,ys,vals,amplitude,10,-1)
	phaseperiod2 = min(diffs2, key=diffs2.get)
	####if outside range search up or down
	n=1
	inc=0.15
	while float(phaseperiod2.split('-')[1]) == (1.0+(n*inc))*period:
		n+=1
		print('initial search > +{0}% range incrementing by +15%'.format(n*inc*100))
		#print('range = {},{}'.format(period*(1+((n-1)*inc)),period*(1+n*inc)))
		diffs2 = serach_sawtooth_range(1+((n-1)*inc),1+n*inc,period,ys,vals,amplitude,10,-1)
		phaseperiod2 = min(diffs2, key=diffs2.get)

	while float(phaseperiod2.split('-')[1]) == (1.0-(n*inc))*period:
		n+=1
		print('initial search < -{0}% range incrementing by -15%'.format(n*inc*100))
		diffs2 = serach_sawtooth_range(1-n*inc,1-((n-1)*inc),period,ys,vals,amplitude,10,-1)
		phaseperiod2 = min(diffs2, key=diffs2.get)
	
	print('initial [phase,period],fitscore',phaseperiod2.split('-'),diffs2[phaseperiod2])
	print('optimizing fit')
	diffs2 = serach_sawtooth_range(0.90,1.10,int(phaseperiod2.split('-')[1]),ys,vals,amplitude,2,-1)
	phaseperiod2 = min(diffs2, key=diffs2.get)
	
	### print top 10 hits for debugging
	dlist = list(diffs2)
	dlist.sort(key=lambda x: diffs2[x])
	for i in dlist[0:11]:
		print('phase: {0:03f} period: {1} score: {2}'.format(float(i.split('-')[0]),int(i.split('-')[1]),diffs2[i]))
	
	if diffs2[phaseperiod2] < diffs[phaseperiod]:
		return(phaseperiod2.split('-'),diffs2[phaseperiod2],-1)
	else:
		return(phaseperiod.split('-'),diffs[phaseperiod],1)






def graph_rot(lines,plotname,maxlength,period):
    HPLs,rots = [],[]
    for i in lines:
        HPLs.append(float(i[labels['_rlnHelicalTrackLength']]))
        rots.append(float(i[labels['_rlnAngleRot']]))
    phaseperiod,difval,sign = eval_sawtooth(180,period,HPLs,rots)
    draw_sawtooth(180,float(phaseperiod[1]),range(0,int(max(HPLs))),float(phaseperiod[0]),sign)
    plt.scatter(HPLs,rots)
    plt.xlim(0,maxlength)
    plt.title('period {0} A, score = {1}'.format(phaseperiod[1],difval))
    plt.savefig('{0}.png'.format(plotname))
    plt.close()
    return([phaseperiod[1],difval])

def graph_xy(lines,plotname):
	xs,ys,mxs,mys = [],[],[],[]
	for i in lines:
		x = float(i[labels['_rlnCoordinateX']])
		y = float(i[labels['_rlnCoordinateY']])
		mx = float(i[labels['_rlnCoordinateX']]) + (float(i[labels['_rlnOriginXAngst']])/1.065)
		my = float(i[labels['_rlnCoordinateY']]) + (float(i[labels['_rlnOriginYAngst']])/1.065)
		plt.scatter(x,y,c='blue')
		plt.scatter(mx,my,c='red')
		plt.plot([x,mx],[y,my],c='black')
	plt.savefig('{0}.png'.format(plotname))
	plt.close()

## read the starfile divide into individual micrographs
errmsg = 'USAGE: rln3p0_HPA.py <starfile> <crossover length>'
try:
	labels,data = read_starfile_new(sys.argv[1])
except:
	sys.exit('\nERROR: Failed to read star file\n{0}'.format(errmsg))
try:
	xover = int(2*(int(sys.argv[2])))
except:
	sys.exit('\nERROR: invalid crossover length\n{0}'.format(errmsg))
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

## iterate over all the mics again for actual calculations
periods,difs = [],[]
finaldatadic = {}				#{micrograph&filament:period,score}
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
		if len(filaments[i]) > 10:
			print('filament:',i,'segments:',len(filaments[i]))
			period,diff = graph_rot(filaments[i],'{0}-{1}_rot'.format(mgname,i),maxlength,xover)
			periods.append(float(period))
			difs.append(float(diff))
			mgfil = '{0}&{1}'.format(mgname,i)
			finaldatadic[mgfil] = [period,diff]
plt.hist(periods,bins=range(int(min(periods)),int(max(periods))+10,10))
print('{0} total fibrils with > 10 segments'.format(len(periods)))
plt.savefig('initial_histogram.png')

output = open('output.txt','w')
for i in finaldatadic:
	mg,fil = i.split('&')
	period,score = (float(finaldatadic[i][0])/2,float(finaldatadic[i][1]))
	output.write('{0},{1},{2},{3}\n'.format(mg,fil,period,score))
output.close()
