from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tree


def main():
	features = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
	with open('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/T9_concatenation.fa','r') as fin:
		ll = fin.readlines()
	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	dic_seq = {i:{} for i in residues}
	for l in ll:
		if l[0] == '>':
			n = l[1:].strip('\n')
		else:
			l = l.strip('\n')
			p = []
			for r in residues:
				dic_seq[r][n] = l.count(r) / (len(l)-l.count('-'))
	df = pd.DataFrame(data = dic_seq)
	
	# Separating out the features
	x = df.loc[:, features].values
	
	pca = PCA()
	principalComponents = pca.fit_transform(x)
	principalDf = pd.DataFrame(data = principalComponents)
	
	
	n = 1
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	ax.set_xlabel('Principal Component 1', fontsize = 15)
	ax.set_ylabel('Principal Component 2', fontsize = 15)
	ax.set_title('2 component PCA', fontsize = 20)
	
	
	
	with open('T9_compositions_4_axes.csv','r') as fin:
		ll = fin.readlines()

	ll = ll[1:]
	dicseq_4 = {}
	for l in ll:
		l = l.strip('\n').split('\t')[:-1]
		dicseq_4[l[0]] = np.array([float(x) for x in l[1:]])


	with open('T9_compositions_0_axes.csv','r') as fin:
		ll = fin.readlines()

	ll = ll[1:]
	dicseq_0 = {}
	for l in ll:
		l = l.strip('\n').split('\t')[:-1]
		dicseq_0[l[0]] = np.array([float(x) for x in l[1:]])
	
	
	Ancestors, BL = tree.main('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/T9.tre')
	des = {}
	for Anc in Ancestors:
		if Ancestors[Anc][0] not in Ancestors and Ancestors[Anc][0] not in des:
			des[Ancestors[Anc][0]] = Anc
		if Ancestors[Anc][1] not in Ancestors and Ancestors[Anc][1] not in des:
			des[Ancestors[Anc][1]] = Anc
		if Ancestors[Anc][2] not in Ancestors and Ancestors[Anc][2] not in des:
			des[Ancestors[Anc][2]] = Anc
	
	coord_4 = {}
	coord_0 = {}
	coord_des = {}
	
# 	max_dist = []
# 	for i,k in enumerate(principalDf.loc[:,:].values):
# 		max_dist.append((k[0]**2 + k[1]**2)**0.5)
# 	max_dist = max(max_dist)
	for i,k in enumerate(principalDf.loc[:,:].values):
# 		ax.scatter([k[0]],[k[1]],c='k')
# 		ax.annotate(df.index[i],(k[0],k[1]))
		coord_des[df.index[i]] = [k[0],k[1]]
	for sp in dicseq_4:
		dicseq_4[sp] = dicseq_4[sp] - pd.DataFrame.mean(df)
# 		ax.scatter([np.sum(pca.components_[0] * dicseq_4[sp])],[np.sum(pca.components_[1] * dicseq_4[sp])],c='r')
		coord_4[sp] = [np.sum(pca.components_[0] * dicseq_4[sp]),np.sum(pca.components_[1] * dicseq_4[sp])]
# 		ax.annotate(sp,(np.sum(pca.components_[0] * dicseq_4[sp]),np.sum(pca.components_[1] * dicseq_4[sp])))
	for sp in dicseq_0:
		dicseq_0[sp] = dicseq_0[sp] - pd.DataFrame.mean(df)
# 		ax.scatter([np.sum(pca.components_[0] * dicseq_0[sp])],[np.sum(pca.components_[1] * dicseq_0[sp])],c='k')
		coord_0[sp] = [np.sum(pca.components_[0] * dicseq_0[sp]),np.sum(pca.components_[1] * dicseq_0[sp])]
# 		ax.annotate(sp,(np.sum(pca.components_[0] * dicseq_0[sp]),np.sum(pca.components_[1] * dicseq_0[sp])))
# 	for i,k in enumerate(zip(pca.components_[0],pca.components_[1])):
# 		ax.plot([0,k[0]*max_dist], [0,k[1]*max_dist])
# 		ax.annotate(features[i],(k[0]*max_dist,k[1]*max_dist))
	
# 	plt.savefig('pca_plots.svg')
	
	for Anc in Ancestors:
		if Anc != 'Anc'+str(len(Ancestors)):
			ax.plot([coord_0[Anc][0],coord_0[Ancestors[Anc][2]][0]],[coord_0[Anc][1],coord_0[Ancestors[Anc][2]][1]],color='k')
			ax.plot([coord_4[Anc][0],coord_4[Ancestors[Anc][2]][0]],[coord_4[Anc][1],coord_4[Ancestors[Anc][2]][1]],color='r')
		ax.scatter([coord_0[Anc][0]],[coord_0[Anc][1]],c='k')
# 		ax.annotate(Anc,(coord_0[Anc][0],coord_0[Anc][1]),c='k')
		ax.scatter([coord_4[Anc][0]],[coord_4[Anc][1]],c='r')
# 		ax.annotate(Anc,(coord_4[Anc][0],coord_4[Anc][1]),c='r')
	for Des in des:
		ax.plot([coord_des[Des][0],coord_0[des[Des]][0]],[coord_des[Des][1],coord_0[des[Des]][1]],color='k')
		ax.plot([coord_des[Des][0],coord_4[des[Des]][0]],[coord_des[Des][1],coord_4[des[Des]][1]],color='r')
		ax.scatter([coord_des[Des][0]],[coord_des[Des][1]],c='k')
		ax.annotate(Des,(coord_des[Des][0],coord_des[Des][1]),c='k')
	plt.show()

	# return pca.components_


if __name__ == "__main__":
	main()