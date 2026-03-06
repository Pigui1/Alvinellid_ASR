import numpy as np
from random import choices


def comp_dic(file_in,dic_ML,Anc_list,dic_tir,tir,th):
	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

	nf = file_in.split('.')[0].replace('_','.')

	with open(file_in,'r') as fin:
		ll = fin.readlines()

	for i,l in enumerate(ll):
		if l.strip('\n') in Anc_list:
			n = nf + '.' + l.strip('\n')
			print(n)
		
			probs = ll[i+1:i+21]
			for j,k in enumerate(probs):
				probs[j] = [float(x) for x in k.strip('\n').split('\t')]
			probs = np.array(probs)
		
			seq_ML = ''
			for k in range(probs.shape[1]):
				seq_ML += residues[list(probs[:,k]).index(max(probs[:,k]))]
			dic_ML[n] = seq_ML
		
			for t in range(tir):
				if (t+1) % (tir/5) <= 1:
					print(t)
				seq = ''
				for j in range(probs.shape[1]):
					if max(probs[:,j]) >= th:
						seq += residues[list(probs[:,j]).index(max(probs[:,j]))]
					else:
						seq += choices(residues,weights=probs[:,j])[0]
				if seq not in dic_tir:
					dic_tir[seq] = {n:1}
				else:
					if n not in dic_tir[seq]:
						dic_tir[seq][n] = 1
					else:
						dic_tir[seq][n] += 1
	return dic_ML,dic_tir



tir = 1000
th = 0.95
dic_ML = {}
dic_tir = {}


file_in = 'MDH_T9_asr_0_axes.csv'
print(file_in)
Anc_list = ['Anc1','Anc2','Anc3','Anc4','Anc5','Anc6','Anc7','Anc8','Anc10','Anc13','Anc16']
dic_ML, dic_tir = comp_dic(file_in,dic_ML,Anc_list,dic_tir,tir,th)


file_in = 'MDH_T9_0missing_model_asr.csv'
print(file_in)
Anc_list = ['Anc1','Anc2','Anc3','Anc4','Anc5','Anc6','Anc7','Anc8','Anc10','Anc13','Anc16']
dic_ML, dic_tir = comp_dic(file_in,dic_ML,Anc_list,dic_tir,tir,th)


file_in = 'MDH_T6_asr_0_axes.csv'
print(file_in)
Anc_list = ['Anc1','Anc2','Anc3','Anc4','Anc5','Anc6','Anc7','Anc9','Anc10','Anc13','Anc16']
dic_ML, dic_tir = comp_dic(file_in,dic_ML,Anc_list,dic_tir,tir,th)

file_in = 'MDH_T6_0missing_model_asr.csv'
print(file_in)
Anc_list = ['Anc1','Anc2','Anc3','Anc4','Anc5','Anc6','Anc7','Anc9','Anc10','Anc13','Anc16']
dic_ML, dic_tir = comp_dic(file_in,dic_ML,Anc_list,dic_tir,tir,th)




with open('ML_seq.fa','w') as fout:
	for n in dic_ML:
		fout.write('>'+n+'\n')
		fout.write(dic_ML[n])
		fout.write('\n')
with open('sampled_seq.fa','w') as fout:
	for s in dic_tir:
		sw = []
		for anc in dic_tir[s]:
			sw.append(anc+'-'+str(dic_tir[s][anc]))
		sw = '_'.join(sw)
		fout.write('>'+sw+'\n')
		fout.write(s)
		fout.write('\n')
