import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
import ML_model
import matrices_model




def get_compositions(compos_mean, coord, eigenvectors):
	mini = 0.0001 #minimum 0.01% to avoid the strict absence of a residue type
	coord = np.array(coord).reshape(-1,1)
	compositions = np.array([x for x in compos_mean]) #initialization with the means
	if coord.size > 0:
		compositions = compositions + np.sum(coord*eigenvectors, axis = 0) #each residue: added coordinate value * residue load on the axis
		if np.min(compositions) <= 0: #correction in a residue % is below 0, and rounding so that sum frequencies == 1
			a = (1-20*mini)/(np.sum(compositions)-20*np.min(compositions))
			b = mini - a*np.min(compositions)
		else:
			a = 1/np.sum(compositions)
			b = 0
		compositions = a*compositions+b
	
	return compositions


def main(df,output):
	cod_aa = 0
	if len(df) == 2:
		cod_aa = 1
		features_nt = ['A','T','C','G','A3','T3','C3','G3']
		df_nt = df[0]
		features_aa = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
		df_aa = df[1]
	elif len(df.loc[:, :].values[0]) == 8:
		features = ['A','T','C','G','A3','T3','C3','G3']
	else:
		features = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
	
	if cod_aa == 0:
		# Separating out the features
		x = df.loc[:, features].values
		
		# Standardizing the features
# 		x = StandardScaler().fit_transform(x)
		
		pca = PCA()
		principalComponents = pca.fit_transform(x)
		principalDf = pd.DataFrame(data = principalComponents)
		
# 		print('explained variance:')
# 		print(100*pca.explained_variance_/np.sum(pca.explained_variance_))
		
		with open(output+'_log_command.txt','a') as f_out:
			f_out.write('explained variance:\n')
			f_out.write('\t'.join([str(x) for x in 100*pca.explained_variance_/np.sum(pca.explained_variance_)])+'\n')
		
	
	if cod_aa == 1:
		x_nt = df_nt.loc[:, features_nt].values
		#x_nt = StandardScaler().fit_transform(x_nt)
		pca_nt = PCA()
		principalComponents_nt = pca_nt.fit_transform(x_nt)
		principalDf_nt = pd.DataFrame(data = principalComponents_nt)
# 		print('explained variance on nucleotides:')
# 		print(100*pca_nt.explained_variance_/np.sum(pca_nt.explained_variance_))
		
		x_aa = df_aa.loc[:, features_aa].values
		#x_aa = StandardScaler().fit_transform(x_aa)
		pca_aa = PCA()
		principalComponents_aa = pca_aa.fit_transform(x_aa)
		principalDf_aa = pd.DataFrame(data = principalComponents_aa)
# 		print('explained variance on amino acids:')
# 		print(100*pca_aa.explained_variance_/np.sum(pca_aa.explained_variance_))
		
		with open(output+'_log_command.txt','a') as f_out:
			f_out.write('explained variance nucleotide:\n')
			f_out.write('\t'.join([str(x) for x in 100*pca_nt.explained_variance_/np.sum(pca_nt.explained_variance_)])+'\n')
			f_out.write('explained variance amino acids:\n')
			f_out.write('\t'.join([str(x) for x in 100*pca_aa.explained_variance_/np.sum(pca_aa.explained_variance_)])+'\n')
	
	
	if len(df) == 2:
		n = 2
	else:
		n = 1
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	ax.set_xlabel('Principal Component 1', fontsize = 15)
	ax.set_ylabel('Principal Component 2', fontsize = 15)
	ax.set_title('2 component PCA', fontsize = 20)
	
	if n == 2:
		max_dist_nt = []
		for i,k in enumerate(principalDf_nt.loc[:,:].values):
			max_dist_nt.append((k[0]**2 + k[1]**2)**0.5)
		max_dist_nt = max(max_dist_nt)
		max_dist_aa = []
		for i,k in enumerate(principalDf_aa.loc[:,:].values):
			max_dist_aa.append((k[0]**2 + k[1]**2)**0.5)
		max_dist_aa = max(max_dist_aa)
		for i,k in enumerate(principalDf_nt.loc[:,:].values):
			ax[0].scatter([k[0]],[k[1]])
			ax[0].annotate(df[0].index[i],(k[0],k[1]))
		for i,k in enumerate(zip(pca_nt.components_[0],pca_nt.components_[1])):
			ax[0].plot([0,k[0]*max_dist_nt], [0,k[1]*max_dist_nt])
			ax[0].annotate(features_nt[i],(k[0]*max_dist_nt,k[1]*max_dist_nt))
		for i,k in enumerate(principalDf_aa.loc[:,:].values):
			ax[1].scatter([k[0]],[k[1]])
			ax[1].annotate(df[1].index[i],(k[0],k[1]))
		for i,k in enumerate(zip(pca_aa.components_[0],pca_aa.components_[1])):
			ax[1].plot([0,k[0]*max_dist_aa], [0,k[1]*max_dist_aa])
			ax[1].annotate(features_aa[i],(k[0]*max_dist_aa,k[1]*max_dist_aa))
	else:
		max_dist = []
		for i,k in enumerate(principalDf.loc[:,:].values):
			max_dist.append((k[0]**2 + k[1]**2)**0.5)
		max_dist = max(max_dist)
		for i,k in enumerate(principalDf.loc[:,:].values):
			ax.scatter([k[0]],[k[1]])
			ax.annotate(df.index[i],(k[0],k[1]))
		for i,k in enumerate(zip(pca.components_[0],pca.components_[1])):
			ax.plot([0,k[0]*max_dist], [0,k[1]*max_dist])
			ax.annotate(features[i],(k[0]*max_dist,k[1]*max_dist))
	
	plt.savefig('pca_plots.svg')
	
	
	if n == 2:
		return pca_nt.components_,pca_aa.components_
	else:
		return pca.components_


def make_plot(df,components_1,components_2,tree,BL,coord_root,model_aa,gam,code,output):
	
	g = ML_model.classesgamma(gam)
	
	if len(df) == 2:
		n = 2
	else:
		n = 1
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	ax.set_xlabel('Principal Component 1', fontsize = 15)
	ax.set_ylabel('Principal Component 2', fontsize = 15)
	ax.set_title('2 component PCA', fontsize = 20)
	
	if n == 2:
		
		nbr_axis = len(components_2)
		
		principalDf_1 = df[0]
		principalDf_2 = df[1]
		
		#no done!!
		
		
		max_dist_1 = []
		for i,k in enumerate(principalDf_1.loc[:,:].values):
			max_dist_1.append((k[0]**2 + k[1]**2)**0.5)
		max_dist_1 = max(max_dist_1)
		max_dist_2 = []
		for i,k in enumerate(principalDf_2.loc[:,:].values):
			max_dist_2.append((k[0]**2 + k[1]**2)**0.5)
		max_dist_2 = max(max_dist_2)
		for i,k in enumerate(principalDf_1.loc[:,:].values):
			ax[0].scatter([k[0]],[k[1]])
			ax[0].annotate(df[0].index[i],(k[0],k[1]))
		for i,k in enumerate(zip(components_1[0],components_1[1])):
			ax[0].plot([0,k[0]*max_dist_1], [0,k[1]*max_dist_1])
			ax[0].annotate(features_1[i],(k[0]*max_dist_1,k[1]*max_dist_1))
		for i,k in enumerate(principalDf_2.loc[:,:].values):
			ax[1].scatter([k[0]],[k[1]])
			ax[1].annotate(df[1].index[i],(k[0],k[1]))
		for i,k in enumerate(zip(components_2[0],components_2[1])):
			ax[1].plot([0,k[0]*max_dist_2], [0,k[1]*max_dist_2])
			ax[1].annotate(features_2[i],(k[0]*max_dist_2,k[1]*max_dist_2))
	else:
		
		nbr_axis = len(components_1)
		
		compo_means = np.array(pd.DataFrame.mean(df))
		comp_init = get_compositions(compo_means, coord_root[1:], components_1)
		compos = {'Anc'+str(len(tree)):comp_init}
				
		for i in range(len(tree),0,-1):
			Anc = 'Anc'+str(i)
			compos[tree[Anc][0]] =  np.sum(np.array([compos[Anc].dot(expm(BL[Anc][0][0] * g[k] * matrices_model.equi_matrix(model_aa,get_compositions(compo_means,BL[Anc][0][1:],components_1),code))) for k in range(4)]),axis=0)/4
			compos[tree[Anc][1]] =  np.sum(np.array([compos[Anc].dot(expm(BL[Anc][1][0] * g[k] * matrices_model.equi_matrix(model_aa,get_compositions(compo_means,BL[Anc][1][1:],components_1),code))) for k in range(4)]),axis=0)/4
			if i == len(tree):
				compos[tree[Anc][2]] =  np.sum(np.array([compos[Anc].dot(expm(BL[Anc][2][0] * g[k] * matrices_model.equi_matrix(model_aa,get_compositions(compo_means,BL[Anc][2][1:],components_1),code))) for k in range(4)]),axis=0)/4
		
		compos_init = {}
		for i,k in enumerate(df.loc[:,:].values):
			compos_init[df.index[i]] = [(k-compo_means).dot(j) for j in components_1]
		compos_simul = {}
		for k in compos:
			compos_simul[k] = [(compos[k]-compo_means).dot(j) for j in components_1]
		
		max_dist = []
		y = 0
		u = 0.01
		for k in compos_init:
			x = compos_init[k][0]
			if nbr_axis > 1:
				y = compos_init[k][1]
				max_dist.append((x**2+y**2)**0.5)
			else:
				y += u
				max_dist.append((x**2)**0.5)
			ax.scatter(x,y,c='black')
			ax.annotate(k,(x,y))
		max_dist = max(max_dist)
		y = 0
		for k in compos_simul:
			x = compos_simul[k][0]
			if nbr_axis > 1:
				y = compos_simul[k][1]
			else:
				y -= u
			ax.scatter(x,y,c='red')
			ax.annotate(k,(x,y))
		if nbr_axis > 1:
			for i,k in enumerate(zip(components_1[0],components_1[1])):
				x = k[0]
				y = k[1]
				ax.plot([0,x*max_dist], [0,y*max_dist], color='black')
				ax.annotate(df.columns[i],(x*max_dist,y*max_dist))
		else:
			for i,k in enumerate(components_1[0]):
				x = k
				y = 0
				ax.plot([0,x*max_dist], [0,y*max_dist], color='black')
				ax.annotate(df.columns[i],(x*max_dist,y*max_dist))	
	
	plt.savefig(output+'_pca_plots_'+str(nbr_axis)+'_axis.svg')
	












