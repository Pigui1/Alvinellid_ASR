import matplotlib.pyplot as plt
import numpy as np
from sys import argv


def get_seq(file_in):
	dic_seq = {}
	with open(file_in,'r') as fin:
		ll = fin.readlines()
	for l in ll:
		if l[0] == '>':
			name = l[1:].strip()
		else:
			dic_seq[name] = l.strip()
	
	length_tot = len(l.strip())
	
	return dic_seq, length_tot


def get_comp(dic_seq, aa_comp, det):
	aa = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	for sp in dic_seq:
		for k in dic_seq[sp]:
			if k in aa:
				aa_comp[sp][aa.index(k)] += 1
			else:
				det[sp] +=1
	
	return aa_comp, det
	

if __name__ == "__main__":
	
	dic_seq, length_tot = get_seq(argv[1])
	
	aa_comp = {x:np.zeros((20,1)) for x in dic_seq}
	det = {x:0 for x in dic_seq}
	aa_comp, det = get_comp(dic_seq, aa_comp, det)
# 	
	for sp in aa_comp:
		aa_comp[sp] /= (length_tot - det[sp])

# 	length_tot1, length_tot2, length_tot3, length_tot4, length_tot5, length_tot6 = 0, 0, 0, 0, 0, 0
# 	
# 	dic_seq1, length_tot1 = get_seq('non_informative_concatenation_1missing_hb.fa')
# 	dic_seq2, length_tot2 = get_seq('non_informative_concatenation_1missing_sb.fa')
# 	dic_seq3, length_tot3 = get_seq('non_informative_concatenation_1missing_tb.fa')
# 	dic_seq4, length_tot4 = get_seq('non_informative_concatenation_1missing_HE.fa')
# 	dic_seq5, length_tot5 = get_seq('non_informative_concatenation_1missing_SE.fa')
# 	dic_seq6, length_tot6 = get_seq('non_informative_concatenation_1missing_TE.fa')
# 	
# 	try:
# 		aa_comp = {x:np.zeros((20,1)) for x in dic_seq1}
# 		det = {x:0 for x in dic_seq1}
# 	except:
# 		aa_comp = {x:np.zeros((20,1)) for x in dic_seq4}
# 		det = {x:0 for x in dic_seq4}
# 	
# 	aa_comp, det = get_comp(dic_seq1, aa_comp, det)
# 	aa_comp, det = get_comp(dic_seq2, aa_comp, det)
# 	aa_comp, det = get_comp(dic_seq3, aa_comp, det)
# 	aa_comp, det = get_comp(dic_seq4, aa_comp, det)
# 	aa_comp, det = get_comp(dic_seq5, aa_comp, det)
# 	aa_comp, det = get_comp(dic_seq6, aa_comp, det)
# 	
# 	for sp in aa_comp:
# 		aa_comp[sp] /= (length_tot1 + length_tot2 + length_tot3 + length_tot4 + length_tot5 + length_tot6 - det[sp])
	
	
	with open('thermostable.txt','r') as fin:
		ll = fin.readlines()
	ll = ll[1:]
	eigenvectors = [[float(y) for y in x.strip().split('\t')] for x in ll]
	
	eigenvectors = np.array(eigenvectors)
	
# 	print(eigenvectors)
	fig = plt.figure(figsize = (8,8))
	ax = fig.add_subplot(1,1,1)
	for sp in aa_comp:
		x = eigenvectors[0].dot(aa_comp[sp])
		y = eigenvectors[1].dot(aa_comp[sp])
		
		ax.scatter(x,y,c='k')
		ax.annotate(sp,(x,y))
	
# 	plt.show()
	plt.savefig(argv[1]+'.svg')
# 	plt.savefig('PCA.svg')




