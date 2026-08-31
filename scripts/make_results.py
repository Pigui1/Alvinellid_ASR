import os
import numpy as np


folder1 = 'long_or_complete'
folder2 = 'others'
folder3 = 'oldgenes'


with open('results_locus.csv','w') as fout:
	fout.write('chromosome\tstart\tstop\tL1\tL2\tL3\tL4\tL5\tL6\tL7\tL8\tL9\tL10\tL11\tL12\tL13\tL14\tL15\tphylogeny\tsites\tlength\torigin\n')


def get_r(folder,r,r_order):
	for root, dirs, files in os.walk('./'+folder+'/', topdown=True):
		for file in files:
			if '.iqtree' in file:
				with open('./'+folder+'/'+file,'r') as fin:
					ll = fin.readlines()

				locus = file.split('-')[0]
				start = file.split('_')[1].split('-')[1]
				stop = file.split('_')[1].split('-')[2]
				phylo = file.split('_')[2].split('.')[0]
				k=0
				for l in ll:
					if l == 'Tree      logL    deltaL\n':
						ind = k
					k+=1
				probs = [x for x in ll[ind+2:ind+17]]
				for i, p in enumerate(probs):
					p = p.strip('\n')
					p = p.split(' ')
					p = [x for x in p if x != '']
					p = p[1]
					probs[i] = p

				sites = 0
				len_align = 0
				alignement = file.split('.')[0]+'.fa'
				with open('./'+folder+'/'+alignement,'r') as fin:
					ll = fin.readlines()
				for l in ll:
					if l[0]!= '>':
						sites += len(l.strip('\n').replace('-',''))
						len_align = len(l.strip('\n'))

			
				if locus not in r:
					r[locus] = []
					r_order[locus] = []
				r[locus].append('\t'.join([locus, start, stop] + probs + [phylo] + [str(sites)] + [str(len_align)] + [folder]+ ['\n']))
				r_order[locus].append(int(start))

		
		return r,r_order



r = {}
r_order = {}
r,r_order = get_r(folder1,r,r_order)
r,r_order = get_r(folder2,r,r_order)
r,r_order = get_r(folder3,r,r_order)


for k in r_order:
	r_order[k] = np.argsort(r_order[k])
	r[k] = [r[k][x] for x in r_order[k]]
	
with open('results_locus.csv','a') as fout:
	for k in r:
		for i in r[k]:
			fout.write(i)







