from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tree


def comp_anc(file_in):
	with open(file_in,'r') as fin:
		ll = fin.readlines()

	ll = ll[1:]
	dicseq = {}
	for l in ll:
		l = l.strip('\n').split('\t')
		dicseq[l[0]] = np.array([float(x) for x in l[1:]])
	
	return dicseq


def main():
	features = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
	with open('./gene_alignments/genes_non_informative/non_informative_concatenation_0missing.fa','r') as fin:
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
	
	# splitting out the features
	x = df.loc[:, features].values
	
	pca = PCA()
	principalComponents = pca.fit_transform(x)
	principalDf = pd.DataFrame(data = principalComponents)
	
	print('explained variance:')
	print(100*pca.explained_variance_/np.sum(pca.explained_variance_))
	
	
	n = 1
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	ax.set_xlabel('Principal Component 1', fontsize = 15)
	ax.set_ylabel('Principal Component 2', fontsize = 15)
	ax.set_title('2 component PCA', fontsize = 20)
	
	
	files_in =  ['T1_0missing_perfect_concat_model_2_axes.csv','T2_0missing_perfect_concat_model_0_axes.csv','T3_0missing_perfect_concat_model_0_axes.csv','T4_0missing_perfect_concat_model_2_axes.csv',
			  'T5_0missing_perfect_concat_model_2_axes.csv','T6_0missing_perfect_concat_model_2_axes.csv','T7_0missing_perfect_concat_model_2_axes.csv','T8_0missing_perfect_concat_model_2_axes.csv',
			  'T9_0missing_perfect_concat_model_2_axes.csv','T10_0missing_perfect_concat_model_2_axes.csv','T11_0missing_perfect_concat_model_2_axes.csv','T12_0missing_perfect_concat_model_2_axes.csv',
			  'T13_0missing_perfect_concat_model_2_axes.csv','T14_0missing_perfect_concat_model_0_axes.csv','T15_0missing_perfect_concat_model_0_axes.csv']
	
# 	files_in = ['T1_0missing_perfect_concat_0axis_model_0_axes.csv', 'T2_0missing_perfect_concat_0axis_model_0_axes.csv', 'T3_0missing_perfect_concat_0axis_model_0_axes.csv', 'T4_0missing_perfect_concat_0axis_model_0_axes.csv', 'T5_0missing_perfect_concat_0axis_model_0_axes.csv',
# 				'T6_0missing_perfect_concat_0axis_model_0_axes.csv', 'T7_0missing_perfect_concat_0axis_model_0_axes.csv', 'T8_0missing_perfect_concat_0axis_model_0_axes.csv', 'T9_0missing_perfect_concat_0axis_model_0_axes.csv', 'T10_0missing_perfect_concat_0axis_model_0_axes.csv',
# 				'T11_0missing_perfect_concat_0axis_model_0_axes.csv','T12_0missing_perfect_concat_0axis_model_0_axes.csv','T13_0missing_perfect_concat_0axis_model_0_axes.csv','T14_0missing_perfect_concat_0axis_model_0_axes.csv','T15_0missing_perfect_concat_0axis_model_0_axes.csv']

	
	
	dicseq_res = [comp_anc(x) for x in files_in]
	
	coord = [{} for x in files_in]
	coord_des = {}
	
	
	for i,k in enumerate(principalDf.loc[:,:].values):
		coord_des[df.index[i]] = [k[0],k[1]]
	
	
	for tree in range(15):
		for sp in dicseq_res[tree]:
			dicseq_res[tree][sp] = dicseq_res[tree][sp] - pd.DataFrame.mean(df)
			coord[tree][sp] = [np.sum(pca.components_[0] * dicseq_res[tree][sp]),np.sum(pca.components_[1] * dicseq_res[tree][sp])]
	
	
	max_dist = []
	for i,k in enumerate(principalDf.loc[:,:].values):
		max_dist.append((k[0]**2 + k[1]**2)**0.5)
	max_dist = max(max_dist)
# 	for i,k in enumerate(principalDf.loc[:,:].values):
# 		ax.scatter([k[0]],[k[1]])
# 		ax.annotate(df.index[i],(k[0],k[1]))
	for i,k in enumerate(zip(pca.components_[0],pca.components_[1])):
		ax.plot([0,k[0]*max_dist], [0,k[1]*max_dist])
		ax.annotate(features[i],(k[0]*max_dist,k[1]*max_dist))
	
	
	
	cl = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']
	labels = ["*","x","d","^","o"]
	
	
	length_align = np.array([8866,6554,7746,11091,14522,21561,20399,25888,26845,12879,15621,14452,14947,7915,7262]) #alignment lengths by topology, from 1 to 15
	mean_pos = [[0,0] for x in range(5)]
	
	#moyenne pondérée par nombre de gènes ou sites
	
	for k in range(15):
		for i,Anc in enumerate([13,14,15,16,17]): #Alvinellidae, Ampharetidae, Terebellidae, Alvi+Ampha, Terebellida ancestors
			Anc_cur = 'Anc' + str(Anc)
			ax.scatter([coord[k][Anc_cur][0]],[coord[k][Anc_cur][1]],c=cl[k],marker=labels[i])
			mean_pos[i][0] += coord[k][Anc_cur][0] * length_align[k] / sum(length_align)
			mean_pos[i][1] += coord[k][Anc_cur][1] * length_align[k] / sum(length_align)
# 			ax.annotate(str(k+1)+'_'+str(Anc),(coord[k][Anc_cur][0],coord[k][Anc_cur][1]),c='k')
	for Des in coord_des:
		ax.scatter([coord_des[Des][0]],[coord_des[Des][1]],c='k')
		ax.annotate(Des,(coord_des[Des][0],coord_des[Des][1]),c='k')
	
	for i,k in enumerate([13,14,15,16,17]):
		ax.scatter([mean_pos[i][0]], [mean_pos[i][1]])
		ax.annotate('mean_'+str(k),(mean_pos[i][0],mean_pos[i][1]))
	
	plt.show()


if __name__ == "__main__":
	main()