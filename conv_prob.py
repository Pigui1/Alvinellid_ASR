import os
import numpy as np
from math import exp


folder1 = 'long_or_complete'
folder2 = 'others'
folder3 = 'oldgenes'


with open('results_locus.csv','r') as fin:
	ll = fin.readlines()

with open('results_locus_probs.csv','w') as fout:
	fout.write(ll[0])
	for l in ll[1:]:
		l = l.split('\t')
		prob = l[3:18]
		prob = [float(x) for x in prob]
		prob = [x-max(prob) for x in prob]
		prob_out = []
		for x in prob:
			if x < -30:
				prob_out.append(0)
			else:
				prob_out.append(exp(x))
		prob_out = np.array(prob_out)/sum(prob_out)
		prob_out = [str(x) for x in prob_out]
		
		l[3:18] = prob_out
		l = '\t'.join(l)
		fout.write(l)