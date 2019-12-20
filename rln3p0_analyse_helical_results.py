#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys

###---------function: read the star file get the header, labels, and data -------------#######

def read_3p1star(infile):
    inDMCloop = False
    labels,data,header = {},[],[]
    dta = open(infile,'r').readlines()
    for i in dta:
        if i.split()== ['data_particles']:
            inDMCloop = True
        elif i.strip('\n') == '# version 30001':
            inDMCloop = False
        if inDMCloop == True:
            if '_rln' in i:
                labels[i.split()[0]] = int(i.split()[1].replace('#',''))-1
                header.append(i)
            elif len(i.split()) == len(labels) and len(i.split()) > 0:
                data.append(i.split())
        if inDMCloop == False:
            header.append(i)
    return(labels,data,header)


#---------------------------------------------------------------------------------------------#

## parse the output file
try:
    data = open(sys.argv[1],'r').readlines()
    labels,stardata,header = read_3p1star(sys.argv[2])
    cutoff = float(sys.argv[3])

except:
    sys.exit('USAGE: analyse_helical_results.py <results file> <orignial star file> <cutoff threshold>')

crossovers = []
for i in data:
    line = i.split(',')
    crossovers.append(float(line[2]))

binz = range(int(min(crossovers)),int(max(crossovers))+5,10)
fig,ax = plt.subplots()
plt.hist(crossovers,bins=binz)
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='down',length=5)
ax.tick_params(axis='x', which='major', direction='down',length=10)

plt.xlabel('crossover length')
plt.ylabel('# filaments')
plt.show()

try:
    cmin = float(raw_input('crossover minimum: '))
    cmax = float(raw_input('crossover maximum: '))
except:
    cmin = float(input('crossover minimum: '))
    cmax = float(input('crossover maximum: '))

output = open('helical_selected.star','w')
for i in header:
    output.write(i)

keepers = []
for i in data:
    line = i.split(',')
    if cmin < float(line[2]) < cmax and float(line[3]) < cutoff:
        keepers.append([line[0],line[1]])
count = 0
for i in stardata:
    id = [i[labels['_rlnMicrographName']].split('/')[-1].split('.')[0],i[labels['_rlnHelicalTubeID']]]
    print(id)
    if id in keepers:
        output.write('{0}\n'.format('   '.join(i)))
        count +=1
print('wrote {0} segemnts from {1} filaments to starfile'.format(count,len(keepers)))
        

