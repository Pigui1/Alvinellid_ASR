from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tree


def coord_nodes(filein,components,means):
	with open(filein,'r') as fin:
		ll = fin.readlines()

	ll = ll[1:]
	dicseq = {}
	for l in ll:
		l = l.strip('\n').split('\t')
# 		length = int(l[-1])
		length = 1
# 		l = l[:-1]
		dicseq[l[0]] = np.array([float(x) for x in l[1:]])	
	
	coord = {}
	
# 	s = np.zeros(20)
# 	for sp in dicseq:
# 		s+=dicseq[sp]
# 	s/=len(dicseq)
	s = means
	
	for sp in dicseq:
		dicseq[sp] = dicseq[sp] - s
		coord[sp] = [np.sum(components[0] * dicseq[sp]),np.sum(components[1] * dicseq[sp])]
		
	Ancestors, _ = tree.main('/Users/pgb/Documents/Data/Roscoff/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/'+filein.split('/')[-1].split('_')[0]+'.tre')
	des = {}
	for Anc in Ancestors:
		if Ancestors[Anc][0] not in Ancestors and Ancestors[Anc][0] not in des:
			des[Ancestors[Anc][0]] = Anc
		if Ancestors[Anc][1] not in Ancestors and Ancestors[Anc][1] not in des:
			des[Ancestors[Anc][1]] = Anc
		if Ancestors[Anc][2] not in Ancestors and Ancestors[Anc][2] not in des:
			des[Ancestors[Anc][2]] = Anc
	
	
	return coord, Ancestors, des, length


def coord_nodes_CvP(filein):
	with open(filein,'r') as fin:
		ll = fin.readlines()

	ll = ll[1:]
	dicseq = {}
	for l in ll:
		l = l.strip('\n').split('\t')
# 		length = int(l[-1])
# 		l = l[:-1]
		length = 1
		dicseq[l[0]] = np.array([float(x) for x in l[1:]])	
	
	coord = {}
	
	
# 	ind1 = 'IVYWREL'
# 	ind2 = 'EDKR'
# 	ind3 = 'GHNPQST'
# 	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
# 	for sp in dicseq:
# # 		dicseq[sp] = dicseq[sp] - s
# 		ind1sp = 0
# 		ind2sp = 0
# 		ind3sp = 0
# 		for r in ind1:
# 			ind1sp += dicseq[sp][residues.index(r)]*length
# 		for r in ind2:
# 			ind2sp += dicseq[sp][residues.index(r)]*length
# 		for r in ind3:
# 			ind3sp += dicseq[sp][residues.index(r)]*length
# 		coord[sp] = [ind1sp,ind2sp-ind3sp]
	
	
	with open('thermostable.txt','r') as fin:
		ll = fin.readlines()
	ll = ll[1:]
	eigenvectors = [[float(y) for y in x.strip().split('\t')] for x in ll]
	eigenvectors = np.array(eigenvectors)
	
	for sp in dicseq:
		x = eigenvectors[0].dot(dicseq[sp])
		y = eigenvectors[1].dot(dicseq[sp])
# 		x = eigenvectors[2].dot(dicseq[sp])
# 		y = eigenvectors[3].dot(dicseq[sp])
		coord[sp] = [x,y]
	
	
	
	Ancestors, _ = tree.main('/Users/pgb/Documents/Data/Roscoff/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/'+filein.split('/')[-1].split('_')[0]+'.tre')
	des = {}
	for Anc in Ancestors:
		if Ancestors[Anc][0] not in Ancestors and Ancestors[Anc][0] not in des:
			des[Ancestors[Anc][0]] = Anc
		if Ancestors[Anc][1] not in Ancestors and Ancestors[Anc][1] not in des:
			des[Ancestors[Anc][1]] = Anc
		if Ancestors[Anc][2] not in Ancestors and Ancestors[Anc][2] not in des:
			des[Ancestors[Anc][2]] = Anc
	
	
	s = np.zeros(2)
	for sp in coord:
		s += coord[sp]
	s/=len(coord)
	
# 	for sp in coord:
# 		coord[sp] = coord[sp]-s
	
	return coord, Ancestors, des, length




def main():
	features = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	cl = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']
	
	with open('/Users/pgb/Documents/Data/Roscoff/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/results_locus_probs.csv','r') as fin:
		ll = fin.readlines()
	ll = ll[1:]
	nbr_genes = np.zeros(15)
	for l in ll:
		l = l.split('\t')[3:18]
		l = [float(x) for x in l]
		if max(l)>0.5:
			nbr_genes[l.index(max(l))]+=1
	
	
	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	dic_seq = {i:{} for i in residues}
	concat = {}
	
# 	folder = './resultat_optim_1missing/'
	folder = ''
	
	
	file_list = ['T1_0missing_burried_perfect_concat_model_2_axes.csv', 'T2_0missing_burried_perfect_concat_model_0_axes.csv', 'T3_0missing_burried_perfect_concat_model_0_axes.csv', 'T4_0missing_burried_perfect_concat_model_2_axes.csv', 'T5_0missing_burried_perfect_concat_model_2_axes.csv',
				 'T6_0missing_burried_perfect_concat_model_2_axes.csv', 'T7_0missing_burried_perfect_concat_model_2_axes.csv', 'T8_0missing_burried_perfect_concat_model_2_axes.csv', 'T9_0missing_burried_perfect_concat_model_2_axes.csv', 'T10_0missing_burried_perfect_concat_model_2_axes.csv',
				 'T11_0missing_burried_perfect_concat_model_2_axes.csv','T12_0missing_burried_perfect_concat_model_2_axes.csv','T13_0missing_burried_perfect_concat_model_2_axes.csv','T14_0missing_burried_perfect_concat_model_0_axes.csv','T15_0missing_burried_perfect_concat_model_0_axes.csv']
	
# 	file_list = ['T1_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T2_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T3_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T4_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T5_0missing_burried_perfect_concat_0axis_model_0_axes.csv',
# 				 'T6_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T7_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T8_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T9_0missing_burried_perfect_concat_0axis_model_0_axes.csv', 'T10_0missing_burried_perfect_concat_0axis_model_0_axes.csv',
# 				 'T11_0missing_burried_perfect_concat_0axis_model_0_axes.csv','T12_0missing_burried_perfect_concat_0axis_model_0_axes.csv','T13_0missing_burried_perfect_concat_0axis_model_0_axes.csv','T14_0missing_burried_perfect_concat_0axis_model_0_axes.csv','T15_0missing_burried_perfect_concat_0axis_model_0_axes.csv']


	
	file_list = [folder+x for x in file_list]
	
# 	k = [int(x.split('_')[0][1:]) for x in file_list]
# 	for i in k:
# 		with open('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/T'+str(i)+'_concatenation.fa','r') as fin:
# 			ll = fin.readlines()
# 		for l in ll:
# 			if l[0] == '>':
# 				n = l[1:].strip('\n')
# 			else:
# 				if n not in concat:
# 					concat[n] = ''
# 				concat[n]+=l.strip('\n')
	
	with open('/Users/pgb/Documents/Data/Roscoff/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes_non_informative/non_informative_concatenation_0missing_burried_residues.fa','r') as fin:
		ll = fin.readlines()
	for l in ll:
		if l[0] == '>':
			n = l[1:].strip('\n')
		else:
			concat[n] = l.strip('\n')
	
	for r in residues:
		for n in concat:
			dic_seq[r][n] = concat[n].count(r) / (len(concat[n])-concat[n].count('-'))
	
	df = pd.DataFrame(data = dic_seq)
	
	# Separating out the features
	x = df.loc[:, features].values
	
	
	pca = PCA()
	principalComponents = pca.fit_transform(x)
	principalDf = pd.DataFrame(data = principalComponents)
# 	print(df)
	print(pca.explained_variance_/np.sum(pca.explained_variance_))
	
	n = 1
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	ax.set_xlabel('Principal Component 1', fontsize = 15)
	ax.set_ylabel('Principal Component 2', fontsize = 15)
	ax.set_title('2 component PCA', fontsize = 20)
	
	
	coord_des = {}
	for i,k in enumerate(principalDf.loc[:,:].values):
		coord_des[df.index[i]] = [k[0],k[1]]
	
	
	max_dist = []
	for i,k in enumerate(principalDf.loc[:,:].values):
		max_dist.append((k[0]**2 + k[1]**2)**0.5)
	max_dist = max(max_dist)
	for i,k in enumerate(zip(pca.components_[0],pca.components_[1])):
		ax.plot([0,k[0]*max_dist], [0,k[1]*max_dist],linewidth=0.8,color='gray')
		ax.annotate(features[i],(k[0]*max_dist,k[1]*max_dist),fontsize=8,color='gray')
	
	
	
	mean_alvi_x = 0
	mean_alvi_y = 0
	mean_ampha_x = 0
	mean_ampha_y = 0
	mean_alvi_ampha_x = 0
	mean_alvi_ampha_y = 0
	mean_tereb_x = 0
	mean_tereb_y = 0
	mean_root_x = 0
	mean_root_y = 0
	len_tot = 0
	size_dot=20
	for i,fl in enumerate(file_list):
		means = pd.DataFrame.mean(df)
		coord,Ancestors,des,length = coord_nodes(fl,pca.components_,means)
		length = nbr_genes[i]
		Anc_Alvi_Ampha = Ancestors[des['Anobothrus']][2]
		if des['Anobothrus'] == Ancestors[Anc_Alvi_Ampha][0]:
			Anc_Alvi = Ancestors[Anc_Alvi_Ampha][1]
		else:
			Anc_Alvi = Ancestors[Anc_Alvi_Ampha][0]
		Anc_Ampha = des['Anobothrus']
		Anc_Tereb = des['Mpalmata']
		Anc_root = des['Pectinaria']
		for Anc in coord:
# 			if Anc != 'Anc17':
# 				ax.plot([coord[Anc][0],coord[Ancestors[Anc][2]][0]],[coord[Anc][1],coord[Ancestors[Anc][2]][1]],color=cl[i])
			if Anc == Anc_Alvi:
				ax.scatter([coord[Anc][0]],[coord[Anc][1]],c=cl[i], s=size_dot, marker='*')
# 				ax.annotate('Alvi_'+str(i+1),(coord[Anc][0],coord[Anc][1]))
				mean_alvi_x += coord[Anc][0]*length
				mean_alvi_y += coord[Anc][1]*length
			if Anc == Anc_Ampha:
				ax.scatter([coord[Anc][0]],[coord[Anc][1]],c=cl[i], s=size_dot, marker='x')
# 				ax.annotate('Ampha_'+str(i+1),(coord[Anc][0],coord[Anc][1]))
				mean_ampha_x += coord[Anc][0]*length
				mean_ampha_y += coord[Anc][1]*length
			if Anc == Anc_Tereb:
				ax.scatter([coord[Anc][0]],[coord[Anc][1]],c=cl[i], s=size_dot, marker='d')
# 				ax.annotate('Tereb_'+str(i+1),(coord[Anc][0],coord[Anc][1]))
				mean_tereb_x += coord[Anc][0]*length
				mean_tereb_y += coord[Anc][1]*length
			if Anc == Anc_root:
				ax.scatter([coord[Anc][0]],[coord[Anc][1]],c=cl[i], s=size_dot, marker='o')
# 				ax.annotate('Root_'+str(i+1),(coord[Anc][0],coord[Anc][1]))
				mean_root_x += coord[Anc][0]*length
				mean_root_y += coord[Anc][1]*length
			if Anc == Anc_Alvi_Ampha:
				ax.scatter([coord[Anc][0]],[coord[Anc][1]],c=cl[i], s=size_dot, marker='^')
# 				ax.annotate('Root_'+str(i+1),(coord[Anc][0],coord[Anc][1]))
				mean_alvi_ampha_x += coord[Anc][0]*length
				mean_alvi_ampha_y += coord[Anc][1]*length
# 		ax.annotate(Anc,(coord_0[Anc][0],coord_0[Anc][1]),c='k')
# 		for Des in des:
# 			ax.plot([coord_des[Des][0],coord[des[Des]][0]],[coord_des[Des][1],coord[des[Des]][1]],color=cl[i])
	
	
	len_tot = np.sum(nbr_genes)
# 	print(nbr_genes)
	mean_alvi_x = mean_alvi_x/(len_tot)
	mean_alvi_y = mean_alvi_y/(len_tot)
	mean_ampha_x = mean_ampha_x/(len_tot)
	mean_ampha_y = mean_ampha_y/(len_tot)
	mean_tereb_x = mean_tereb_x/(len_tot)
	mean_tereb_y = mean_tereb_y/(len_tot)
	mean_root_x = mean_root_x/(len_tot)
	mean_root_y = mean_root_y/(len_tot)
	mean_alvi_ampha_x = mean_alvi_ampha_x/(len_tot)
	mean_alvi_ampha_y = mean_alvi_ampha_y/(len_tot)
	
	ax.scatter([mean_alvi_x],[mean_alvi_y],c='k')
	ax.annotate('MEAN_ALVI',(mean_alvi_x,mean_alvi_y))
	ax.scatter([mean_ampha_x],[mean_ampha_y],c='k')
	ax.annotate('MEAN_AMPHA',(mean_ampha_x,mean_ampha_y))
	ax.scatter([mean_tereb_x],[mean_tereb_y],c='k')
	ax.annotate('MEAN_TEREB',(mean_tereb_x,mean_tereb_y))
	ax.scatter([mean_root_x],[mean_root_y],c='k')
	ax.annotate('MEAN_ROOT',(mean_root_x,mean_root_y))
	ax.scatter([mean_alvi_ampha_x],[mean_alvi_ampha_y],c='k')
	ax.annotate('MEAN_ALVI_AMPHA',(mean_alvi_ampha_x,mean_alvi_ampha_y))
	
	
	for Des in des:
		ax.scatter([coord_des[Des][0]],[coord_des[Des][1]],c='k')
		ax.annotate(Des,(coord_des[Des][0],coord_des[Des][1]),c='k')
	
	
# 	plt.show()
# 	exit()
	
	plt.clf()
	
	
	
	
	fig = plt.figure(figsize = (n*8,8))
	ax = fig.add_subplot(1,n,1)
	
	
	coord_tot_des = {}
# 	ind1 = 'IVYWREL'
# 	ind2 = 'EDKR'
# 	ind3 = 'GHNPQST'
	
# 	for n in concat:
# 		coord_tot_des[n] = [0,0]
# 		for r in ind1:
# 			coord_tot_des[n][0] += concat[n].count(r)
# 		coord_tot_des[n][0] = coord_tot_des[n][0] / (len(concat[n])-concat[n].count('-'))
# 		for r in ind2:
# 			coord_tot_des[n][1] += concat[n].count(r)
# 		for r in ind3:
# 			coord_tot_des[n][1] -= concat[n].count(r)
# 		coord_tot_des[n][1] = coord_tot_des[n][1] / (len(concat[n])-concat[n].count('-'))
	
	
	
	
	with open('thermostable.txt','r') as fin:
		ll = fin.readlines()
	ll = ll[1:]
	eigenvectors = [[float(y) for y in x.strip().split('\t')] for x in ll]
	eigenvectors = np.array(eigenvectors)
	
	for sp in concat:
		sp_count = np.array([concat[sp].count(r) for r in residues]) / (len(concat[sp])-concat[sp].count('-'))
		x = eigenvectors[0].dot(sp_count)
		y = eigenvectors[1].dot(sp_count)
# 		x = eigenvectors[2].dot(sp_count)
# 		y = eigenvectors[3].dot(sp_count)
		coord_tot_des[sp] = [x,y]
	
	
	
	s = np.zeros(2)
	for sp in coord_tot_des:
		s += coord_tot_des[sp]
	s/=len(coord_tot_des)
	
# 	for sp in coord_tot_des:
# 		coord_tot_des[sp] = coord_tot_des[sp]-s
	
	len_tot = 0
	coord_tot = {'Anc_uni':[],'Anc_alvi':[],'Anc_pando':[],'Anc_para2':[],'Anc_pomp2':[]}
	for fl in file_list:
		i = int(fl.split('/')[-1].split('_')[0][1:])-1
# 		means = pd.DataFrame.mean(df)
		coord,Ancestors,des,length = coord_nodes_CvP(fl)
		length = nbr_genes[i]
		Anc_conv = {}
		Anc_conv[des['Pectinaria']] = 'Anc_tot'
		Anc_conv[des['Mpalmata']] = 'Anc_palm'
		Anc_conv[des['Skaia']] = 'Anc_kaia'
		Anc_conv[des['Anobothrus']] = 'Anc_ano'
		Anc_conv[Ancestors[des['Anobothrus']][2]] = 'Anc_ampha_alvi'
		Anc_conv[des['Hinvalida']] = 'Anc_inva'
		Anc_conv[des['Acarldarei']] = 'Anc_carlda'
		Anc_conv[des['Pnov']] = 'Anc_nov'
		Anc_conv[des['Psulfincola']] = 'Anc_sulf'
		Anc_conv[des['Phessleri']] = 'Anc_hess'
		Anc_conv[Ancestors[des['Phessleri']][2]] = 'Anc_hess2'
		Anc_conv[des['Pgrasslei']] = 'Anc_grass'
		Anc_conv[Ancestors[des['Pgrasslei']][2]] = 'Anc_para'
		Anc_conv[des['Apompejana']] = 'Anc_pomp'
		if des['Anobothrus'] == Ancestors[Ancestors[des['Anobothrus']][2]][0]:
			k = Ancestors[Ancestors[des['Anobothrus']][2]][1]
			Anc_alvi = k
# 			Anc_conv[k] = 'Anc_alvi'
			Anc_conv[k] = 0
		else:
			k = Ancestors[Ancestors[des['Anobothrus']][2]][0]
			Anc_alvi = k
			Anc_conv[k] = 0
		Anc_uni = des['Punidentata']
		Anc_conv[des['Punidentata']] = 0
		Anc_pando = des['Ppandorae']
		Anc_conv[des['Ppandorae']] = 0
		k = Ancestors[des['Pgrasslei']][2]
		Anc_para2 = Ancestors[k][2]
		Anc_conv[Ancestors[k][2]] = 0
		Anc_pomp2 = Ancestors[des['Apompejana']][2]
		Anc_conv[Ancestors[des['Apompejana']][2]] = 0
# 		len_tot += length
		for Anc in coord:
			if Anc_conv[Anc] not in coord_tot and Anc_conv[Anc] != 0:
				coord_tot[Anc_conv[Anc]] = [0,0]
			if Anc == Anc_alvi:
				coord_tot['Anc_alvi'].append([x for x in coord[Anc]])
			if Anc == Anc_uni:
				coord_tot['Anc_uni'].append([x for x in coord[Anc]])
			if Anc == Anc_pando:
				coord_tot['Anc_pando'].append([x for x in coord[Anc]])
			if Anc == Anc_para2:
				coord_tot['Anc_para2'].append([x for x in coord[Anc]])
			if Anc == Anc_pomp2:
				coord_tot['Anc_pomp2'].append([x for x in coord[Anc]])
			if Anc != Anc_alvi and Anc != Anc_uni and Anc != Anc_pando and Anc != Anc_para2 and Anc != Anc_pomp2:
				coord_tot[Anc_conv[Anc]][0] += coord[Anc][0]*length
				coord_tot[Anc_conv[Anc]][1] += coord[Anc][1]*length
	
	len_tot = np.sum(nbr_genes)
	for Anc in coord_tot:
		if Anc != 'Anc_alvi' and Anc != 'Anc_uni' and Anc != 'Anc_pando' and Anc != 'Anc_para2' and Anc != 'Anc_pomp2':
			coord_tot[Anc][0] = coord_tot[Anc][0]/len_tot
			coord_tot[Anc][1] = coord_tot[Anc][1]/len_tot
	
	c = 'k'
	x = coord_tot_des['Pectinaria']
	y = coord_tot['Anc_tot']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Skaia']
	y = coord_tot['Anc_kaia']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Nedwardsi']
	y = coord_tot['Anc_kaia']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Mpalmata']
	y = coord_tot['Anc_palm']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Acarldarei']
	y = coord_tot['Anc_carlda']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Agunneri']
	y = coord_tot['Anc_inva']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Hinvalida']
	y = coord_tot['Anc_inva']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Anobothrus']
	y = coord_tot['Anc_ano']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Pfijiensis']
	y = coord_tot['Anc_nov']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Pnov']
	y = coord_tot['Anc_nov']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Psulfincola']
	y = coord_tot['Anc_sulf']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Phessleri']
	y = coord_tot['Anc_hess']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Pmira']
	y = coord_tot['Anc_hess']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Ppalmiformis']
	y = coord_tot['Anc_grass']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Pgrasslei']
	y = coord_tot['Anc_grass']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Apompejana']
	y = coord_tot['Anc_pomp']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot_des['Acaudata']
	y = coord_tot['Anc_pomp']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_kaia']
	y = coord_tot['Anc_palm']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_palm']
	y = coord_tot['Anc_tot']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_inva']
	y = coord_tot['Anc_carlda']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_carlda']
	y = coord_tot['Anc_ano']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_ano']
	y = coord_tot['Anc_ampha_alvi']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_ampha_alvi']
	y = coord_tot['Anc_tot']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_nov']
	y = coord_tot['Anc_sulf']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_sulf']
	y = coord_tot['Anc_hess2']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_hess']
	y = coord_tot['Anc_hess2']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_hess2']
	y = coord_tot['Anc_para']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = coord_tot['Anc_grass']
	y = coord_tot['Anc_para']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	for i in range(15):
		x = coord_tot['Anc_alvi'][i]
		plt.scatter([x[0]],[x[1]],c=cl[i],marker='d')
		x = coord_tot['Anc_uni'][i]
		plt.scatter([x[0]],[x[1]],c=cl[i],marker='^')
		x = coord_tot['Anc_pando'][i]
		plt.scatter([x[0]],[x[1]],c=cl[i],marker='o')
		x = coord_tot['Anc_para2'][i]
		plt.scatter([x[0]],[x[1]],c=cl[i],marker='x')
		x = coord_tot['Anc_pomp2'][i]
		plt.scatter([x[0]],[x[1]],c=cl[i],marker='+')
	x_alvi = np.array([0.0,0.0])
	x_uni = np.array([0.0,0.0])
	x_pando = np.array([0.0,0.0])
	x_para2 = np.array([0.0,0.0])
	x_pomp2 = np.array([0.0,0.0])
	for i in range(15):
		x_alvi += (np.array(coord_tot['Anc_alvi'][i]) * nbr_genes[i] / np.sum(nbr_genes))
		x_pando += (np.array(coord_tot['Anc_pando'][i]) * nbr_genes[i] / np.sum(nbr_genes))
		x_uni += (np.array(coord_tot['Anc_uni'][i]) * nbr_genes[i] / np.sum(nbr_genes))
		x_para2 += (np.array(coord_tot['Anc_para2'][i]) * nbr_genes[i] / np.sum(nbr_genes))
		x_pomp2 +=( np.array(coord_tot['Anc_pomp2'][i]) * nbr_genes[i] / np.sum(nbr_genes))
	x = x_alvi
	y = coord_tot['Anc_ampha_alvi']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = x_pando
	y = coord_tot_des['Ppandorae']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = x_uni
	y = coord_tot_des['Punidentata']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = x_para2
	y = coord_tot['Anc_para']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)
	x = x_pomp2
	y = coord_tot['Anc_pomp']
	plt.plot([x[0],y[0]],[x[1],y[1]],color=c)	
	
		
		
# 		y = coord_tot['Anc_ampha_alvi']
# 		plt.scatter([x[0],y[0]],[x[1],y[1]],c=cl[i])
	
	
	
# 	for Anc in coord_tot:
# 		if Anc != 'Anc_tot':
# 			ax.plot([coord_tot[Anc][0],coord_tot[Ancestors[Anc][2]][0]],[coord_tot[Anc][1],coord_tot[Ancestors[Anc][2]][1]],color='k')
	for Des in des:
# 		ax.plot([coord_tot_des[Des][0],coord_tot[des[Des]][0]],[coord_tot_des[Des][1],coord_tot[des[Des]][1]],color=cl[i])
		ax.annotate(Des,(coord_tot_des[Des][0],coord_tot_des[Des][1]))
	plt.show()




#test especes contemporaines	
# 		with open('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/analyse_proteome/gene_alignments/genes/'+fl.split('_')[0]+'_concatenation.fa','r') as fin:
# 			ll = fin.readlines()
# 		residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
# 		prob = {}
# 		for l in ll:
# 			if l[0] == '>':
# 				n = l[1:].strip('\n')
# 			else:
# 				l = l.strip('\n')
# 				p = []
# 				for r in residues:
# 					p.append(l.count(r) / (len(l)-l.count('-')))
# 				prob[n] = np.array([x for x in p])
# 		s = np.zeros(20)
# 		for n in prob:
# 			s+=prob[n]
# 		s/=20
# 		for n in prob:
# 			p = prob[n]-s
# # 			p = np.array(p) - means
# 			c_temp = [np.sum(pca.components_[0] * p),np.sum(pca.components_[1] * p)]
# 			if n =='Apompejana' or n == 'Ppandorae':
# 				ax.scatter([c_temp[0]],[c_temp[1]],c=cl[i])
# 				ax.annotate(n,(c_temp[0],c_temp[1]))
	
	

	
# 	plt.savefig('pca_plots.svg')


	


if __name__ == "__main__":
	main()