import os
import numpy as np
from math import exp, isfinite
import matplotlib.pyplot as plt
from random import sample
from copy import deepcopy
import pandas as pd
from subprocess import run, PIPE
from scipy.stats import binom
import multiprocessing as mp
from random import random


def get_duplicates_stats(r,ws,dis):
	ind = np.array(sample(range(len(r)),k=len(r)))
	
	pval_out = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] #p-values obtenues sur les fenêtres pour toutes les topologies
	
	for chromo in range(1,18):
		max_pos = max(r['stop'].loc[chromo]) - ws/2
		pos = min(ws/2, max_pos)
		
		mv = ws/10 #une fenêtre calculée tous les xxx nucléotides
		
		while pos <= max_pos:
			
			mask = (r['start'] >= pos-ws/2) & (r['stop'] <= pos+ws/2) & (r.index==chromo)
			ss = pd.DataFrame.sum(r.iloc[ind[mask]].iloc[:,22:37], axis=0)
			
			proba = np.array(ss) #probabilités strictes
			
			if len(r.loc[mask]) > 0:
				b = binom(len(r.loc[mask]),dis)
				
				sf = b.sf(proba)
				pvals = np.where(sf < 0.5, b.logsf(proba), -b.logcdf(proba))
				
				for i,p in enumerate(pvals):
					pval_out[i].append(p)
			
			pos += mv
			if pos < max_pos+mv:
				pos = min(pos,max_pos)
	
	for i in range(15):
		finite_vals = [x for x in pval_out[i] if isfinite(x)]
		p_max = max(finite_vals)
		p_min = min(finite_vals)
		pval_out[i] = [p_max if x == float("inf") else p_min if x == float("-inf") else x for x in pval_out[i]] #clipping values
	
	return pval_out
			


def get_res(r,ws,dis):
	
	ret = {}
	pval = {}
	
	for chromo in range(1,18):
		ret[chromo] = []
		pval[chromo] = []
		phylo = -1
		max_pos = max(r['stop'].loc[chromo]) - ws/2
		pos = min(ws/2, max_pos)
		
		mv = ws/10 #une fenêtre calculée tous les xxx nucléotides
		
		while pos <= max_pos:
		
			mask = (r['start'] >= pos-ws/2) & (r['stop'] <= pos+ws/2) & (r.index==chromo)
			ss = pd.DataFrame.sum(r.loc[mask].iloc[:,22:37], axis=0)
			l_ali = pd.DataFrame.sum(r.loc[mask,'length'], axis=0)
			
			proba = np.array(ss) #probabilités strictes.
			
			if len(ret[chromo]) == 0:
				ret[chromo].append([1,list(proba).index(max(proba)),l_ali]) #pseudo-count for display only
			
			ret[chromo].append([pos,list(proba).index(max(proba)),l_ali])

			if len(r.loc[mask]) > 0:

				b = binom(len(r.loc[mask]),dis)
				pw = b.logsf(ss)
				
				sf = b.sf(proba)
				pw = np.where(sf < 0.5, b.logsf(proba), -b.logcdf(proba))
				
				pval[chromo].append([pos,pw])
			pos += mv
			if pos < max_pos+mv:
				pos = min(pos,max_pos)
		
		ret[chromo].append([max_pos + ws/2, list(proba).index(max(proba)),l_ali]) #pseudo-count for display only
	
	return ret,pval



def get_res_3topos(r,ws,dis):
	
	ret = {}
	pval = {}
	
	for chromo in range(1,18):
		ret[chromo] = []
		max_pos = max(r['stop'].loc[chromo]) - ws/2
		pos = min(ws/2, max_pos)
		
		mv = ws/10 #une fenêtre calculée tous les xxx nucléotides
		
		while pos <= max_pos:
		
			mask = (r['start'] >= pos-ws/2) & (r['stop'] <= pos+ws/2) & (r.index==chromo)
			ss = pd.DataFrame.sum(r.loc[mask].iloc[:,2:17], axis=0)
			l_ali = pd.DataFrame.sum(r.loc[mask,'length'], axis=0)
			
			proba = np.array(ss) #probabilités strictes.
			topo6 = proba[2] + proba[5] + proba[7] + proba[10] + proba[12]
			topo7 = proba[1] + proba[3] + proba[6] + proba[11] + proba[14]
			topo9 = proba[0] + proba[4] + proba[8] + proba[9] + proba[13]
			
			proba[2] = topo6
			proba[5] = topo6
			proba[7] = topo6
			proba[10] = topo6
			proba[12] = topo6
			proba[1] = topo7
			proba[3] = topo7
			proba[6] = topo7
			proba[11] = topo7
			proba[14] = topo7
			proba[0] = topo9
			proba[4] = topo9
			proba[8] = topo9
			proba[9] = topo9
			proba[13] = topo9
			
			
			if len(ret[chromo]) == 0:
				ret[chromo].append([1,list(proba).index(max(proba)),l_ali]) #pseudo-count for display only
			
			ret[chromo].append([pos,list(proba).index(max(proba)),l_ali])
			
			pos += mv
			if pos < max_pos+mv:
				pos = min(pos,max_pos)
		
		ret[chromo].append([max_pos + ws/2, list(proba).index(max(proba)),l_ali]) #pseudo-count for display only
	
	return ret	
	

def get_prob(v):
	v2 = v-np.max(v)
	v2 = np.array([max(-30,x) for x in v2])
	v2 = np.exp(v2)
	return v2/np.sum(v2)


if __name__ == '__main__':
	p = 0.0 #filtre sur la probabilité minimale des gènes. 0.33, 0.5, 0.8
	pr_al = 0.0 #filtre sur la proportion minimum de sites alignés par gène
	ws = 2.5e6 #taille de la fenêtre sur le génome
# 	ws = 3.5e7 #taille alternative de la fenêtre sur le génome
	dup0 = 1000 #nombre de réplicats aléatoires, minimum 20 pour p-value = 0.05, 100 pour p-value = 0.01
# 	dup0 = 100


	r = []
	nbr_genes = 0
	with open('results_locus_probs.csv','r') as fin:
		l = fin.readline()
		col = l.strip('\n').split('\t')[1:]+['prop_ali']+['T'+str(i) for i in range(1,16)]
		ind = []
		for counter, l in enumerate(fin.readlines()):
			l = l.strip('\n').split('\t')
			l[0] = l[0].split('_')[0]
			l[1] = int(l[1])
			l[2] = int(l[2])
			l[3:18] = [float(x) for x in l[3:18]]
			m = max(l[3:18])
			c = l[3:18].count(m)
			lint = []
			for k in l[3:18]:
				if k == m and k > p:
					lint.append(1/c)
				else:
					lint.append(0)
			l[19] = int(l[19])
			l[20] = int(l[20])
			if l[18] == 'A':
				nbr_br = 9
			elif l[18] in 'FCB':
				nbr_br = 8
			elif l[18] in 'GHD':
				nbr_br = 7
			elif l[18] in 'EI':
				nbr_br = 6
			elif l[18] in 'KJ':
				nbr_br = 5
			elif l[18] in 'L':
				nbr_br = 4
			pg = l[19]/(nbr_br*l[20])
			if sum(lint)>0:
				ind.append(int(l[0][7:11]))
				r.append(l[1:]+[pg]+lint)
				nbr_genes += 1
	print(nbr_genes)

	r = pd.DataFrame(r,columns=col,index=ind)
	r = r.loc[r['prop_ali']>pr_al]

	col = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']


#################################################################################
##### Figure % of topology given filter threshold, with confidence intervals ####
#################################################################################
# 	p = 0.05
# 	q=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 	q_conf_1=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 	q_conf_2=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 	qtot=[]
# 	x=[]
# 	maxprop = 0
# 	
# 	while p<=0.95:
# 		qinter = []
# 		for k in range(1,16):
# 			L = 'L'+str(k)
# 			T = 'T'+str(k)
# 			s = r.loc[r[L]>p] #should be above the confidence threshold
# # 			s = s.loc[s[L] < p+0.25]
# 			s = s.loc[s[T] > 0] #should be the maximum likelihood topology
# 			qinter.append(np.sum(s[T])) #number of genes under that topology
# 		for k in range(15):
# 			prop_inter = qinter[k]/sum(qinter)
# 			n = sum(qinter)
# 			st = 1.96 * ((0.5*1.96/n)**2 + prop_inter*(1-prop_inter) / n)**0.5
# 			conf_1 = (prop_inter + 1.96**2 / (2*n) + st) / (1 + (1.96**2)/n)
# 			conf_2 = (prop_inter + 1.96**2 / (2*n) - st) / (1 + (1.96**2)/n)
# 			
# 			q[k].append(prop_inter)
# 			q_conf_1[k].append(conf_1)
# 			q_conf_2[k].append(conf_2)
# 		qtot.append(sum(qinter))
# 		x.append(p)
# 		p+=0.01
# 	
# 	fig, ax = plt.subplots()
# 
# 	for i, (qi, qi1, qi2) in enumerate(zip(q, q_conf_1, q_conf_2)):
# 		if i in [5, 8]:
# 			ax.fill_between(x, qi1, qi2, color=col[i], alpha=0.3)
# 
# 	maxprop = 0
# 	for i, qi in enumerate(q):
# 		ax.plot(x, qi, color=col[i])
# 		ax.annotate(f"T{i+1}", (x[-1], qi[-1]), color=col[i])
# 		maxprop = max(max(qi), maxprop)
# 
# 	scale = maxprop / qtot[0]
# 	ax.plot(x, np.array(qtot) * scale, color='k')
# 
# 	ax2 = ax.twinx()
# 	ax2.set_ylim(np.array(ax.get_ylim()) / scale)
# 	
# 	plt.show()
# 	
# 	exit()
#################################################################################
##### Figure % of topology given filter threshold, with confidence intervals ####
#################################################################################




###########################################################################
##### distribution of probabilities, needed for p_values_histograms.py ####
###########################################################################
# 	dis = []
# 	for k in range(1,16):
# 		L = 'L'+str(k)
# 		T = 'T'+str(k)
# 		s = r.loc[r[L]>p]
# 		s = s.loc[s[T] > 0]
# #  		dis.append(len(s))
# 		dis.append(np.sum(s[T]))
# 	dis = [x/sum(dis) for x in dis]
# 	print(dis)
# 	
# 	ret,pval = get_res(r,ws,dis)
# 	pval_tot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 	
# 	for chr in pval:
# 		for loc in pval[chr]:
# 			for i,_ in enumerate(pval_tot):
# 				pval_tot[i].append(loc[1][i])
# 	for i,_ in enumerate(pval_tot):
# 		pval_tot[i] = sorted(pval_tot[i])
# 		print(i,pval_tot[i])
# 	exit()
###########################################################################
##### distribution of probabilities, needed for p_values_histograms.py ####
###########################################################################	
	

	dis = []
	for k in range(1,16):
		L = 'L'+str(k)
		T = 'T'+str(k)
		s = r.loc[r[L]>p]
		s = s.loc[s[T] > 0]
		dis.append(np.sum(s[T]))
	dis = [x/sum(dis) for x in dis]
	print(dis)

	pval_tot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
	pval_tot_util = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
	dup = int(round(ws/1e6 * dup0)) #adapter le nombre de simulations en fonction du nombre de fenêtres dans le génome
	
	
	
	with mp.Pool(processes=4) as pool:
		pvals_out = pool.starmap(get_duplicates_stats,[(r,ws,dis)]*dup)
	for pval in pvals_out:
		for i,_ in enumerate(pval_tot):
			pval_tot[i]+=pval[i]
			pval_tot_util[i].append(min(pval[i])) #meilleure p-value obtenue par une fenêtre générée aléatoirement, où la concentration d'une topologie est la plus élevée
	
	
	pvals = [-1]*15
	for k,p in enumerate(pval_tot_util):
		p = sorted(p)
		pvals[k] = p[int(dup*0.0242-1)] #0.0242 : p-value 5% for 15 independant tested totpologies
	
	
	pvals1 = [-1]*15
	pvals10 = [-1]*15
	for k,p in enumerate(pval_tot_util):
		p = sorted(p)
		pvals1[k] = p[int(dup*0.01-1)] #0.01 : p-value 1% for 15 independant tested totpologies
		pvals10[k] = p[int(dup*0.036-1)] #0.036 : p-value 10% for 15 independant tested totpologies
		
	
	lambda_pvals = np.array([-1]*15)
	for k,p in enumerate(pval_tot_util):
		p = np.exp(p)
		lambda_pvals[k] = 1 / np.mean(p) #can be fitted to an exponential law
	
	print(dup)
	print(pvals)
	print(lambda_pvals)
	
	with open('../figures_resultats/result_simulations_'+str(ws)+'.csv','w') as fout:
		for i in range(15):
			fout.write('\t'.join([str(x) for x in sorted(pval_tot_util[i])]) + '\n')
	
	
	ret,pval = get_res(r,ws,dis)
	ret3 = get_res_3topos(r,ws,dis)
	
	
	
	###Code for global imbalence test
	
	L_exp_global = []
	min_th = -np.log(0.5)
	for t in range(15):
		L_inter = []
		for ch in pval:
			for w in pval[ch]:
				L_inter.append(np.abs(w[1][t]))
		mean_L_inter = sum(L_inter) / len(L_inter) - min_th
		if not isfinite(mean_L_inter):
			mean_L_inter = 10.0
		L_exp_global.append(mean_L_inter)
	
	
	nb_window = len(pvals_out[0][0])
	pvals_windows = [[[] for _ in range(nb_window)] for _ in range(15)] #topology, nb window, replicate
	for k in range(dup):
		for i,_ in enumerate(pvals_windows): #each topology i
			for j in range(nb_window): #each window
				pvals_windows[i][j].append(pvals_out[k][i][j])
	
	fig, axes = plt.subplots(
	nrows=5,
	ncols=3,
	figsize=(15, 12),
	sharex=True,
	sharey=True
	)
	axes = axes.ravel()
	
	cl = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']
	for ax, t in zip(axes,range(15)):
		L_genome_sim = [[0 for _ in range(nb_window)] for _ in range(dup)]
		for u,k in enumerate(pvals_windows[t]): #windows
			for i,j in enumerate(k): #replicate
				L_genome_sim[i][u] = np.abs(j)
		
		for i,_ in enumerate(L_genome_sim):
			L_genome_sim[i] = np.mean(L_genome_sim[i]) - min_th
		ax.hist(L_genome_sim, bins=int(dup**0.5),color=cl[t])
		ax.plot([L_exp_global[t],L_exp_global[t]],[0,len(L_genome_sim)/5],color='k',lw=1)
		
		counter = 0
		for l2 in L_genome_sim:
			if l2 > L_exp_global[t]:
				counter += 1
		pval_genome = counter / len(L_genome_sim)
		ax.set_title(f"Topology {t+1}, pval: {pval_genome}")
	
	
	plt.savefig('../figures_resultats/chromosomes_imbalance_genome'+str(ws)+'.svg') # topologies with imbalanced localization across the genome
	

	cl_3topos_dic = {0:8,1:6,2:5,3:6,4:8,5:5,6:6,7:5,8:8,9:8,10:5,11:6,12:5,13:8,14:6}
	
	fig, axs = plt.subplots(ncols=5, nrows=4, sharey=True, sharex=True, figsize=(3.75, 10), layout="constrained")
	lw = 0.5
	d_max = 0
	chrom_max = 0
	for k in range(1,18):
		for j in ret[k]:
			d_max = max(d_max,j[2])
			if d_max == j[2]:
				chrom_max = k
	print(chrom_max,d_max)
	k = 1
	for row in range(4):
		for col in range(5):
			if k <=17:
				y_n = 1
				y_d = [0]
				x_d = [0]
				for i,j in enumerate(ret[k][:-1]):
					y_d.append(j[0])
					x_d.append(-j[2]/d_max)
					y1 = j[0]
					y2 = (j[0]+ret[k][i+1][0])/2
					y3 = ret[k][i+1][0]
					c1 = cl[j[1]]
					c2 = cl[ret[k][i+1][1]]
					if c1 != c2:
						axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [y2,y2],color=c1)
						y_n = y2
						axs[row,col].plot([0.0,1.0],[y2,y2],color='k',linewidth=lw)
				axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [ret[k][-1][0],ret[k][-1][0]],color=c1)
				axs[row,col].plot([0.0,1.0],[1,1],color='k',linewidth=lw)
				axs[row,col].plot([0.0,1.0],[ret[k][-1][0],ret[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([0.0,0.0],[1,ret[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([1.0,1.0],[1,ret[k][-1][0]],color='k',linewidth=lw)
				x_d.append(-ret[k][-1][2]/d_max)
				y_d.append(ret[k][-1][0])
				x_d.append(0)
				y_d.append(ret[k][-1][0])
				axs[row,col].plot(x_d,y_d,'k',linewidth=lw)
				
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
				
				k+=1
	
	plt.savefig('../figures_resultats/chromosomes_'+str(ws)+'.svg') #coloring according to 15 topologies original probabilities

	fig, axs = plt.subplots(ncols=5, nrows=4, sharey=True, sharex=True, figsize=(3.75, 10), layout="constrained")
	k = 1
	for row in range(4):
		for col in range(5):
			if k <=17:
				y_n = 1
				y_d = [0]
				x_d = [0]
				for i,j in enumerate(ret[k][:-1]):
					y_d.append(j[0])
					x_d.append(-j[2]/d_max)
					y1 = j[0]
					y2 = (j[0]+ret[k][i+1][0])/2
					y3 = ret[k][i+1][0]
					c1 = cl[cl_3topos_dic[j[1]]] #collapse on 3 topos
					c2 = cl[cl_3topos_dic[ret[k][i+1][1]]] #collapse on 3 topos
					if c1 != c2:
						axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [y2,y2],color=c1)
						y_n = y2
						
						axs[row,col].plot([0.0,1.0],[y2,y2],color='k',linewidth=lw)
				axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [ret[k][-1][0],ret[k][-1][0]],color=c1)
				axs[row,col].plot([0.0,1.0],[1,1],color='k',linewidth=lw)
				axs[row,col].plot([0.0,1.0],[ret[k][-1][0],ret[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([0.0,0.0],[1,ret[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([1.0,1.0],[1,ret[k][-1][0]],color='k',linewidth=lw)
				x_d.append(-ret[k][-1][2]/d_max)
				y_d.append(ret[k][-1][0])
				x_d.append(0)
				y_d.append(ret[k][-1][0])
				axs[row,col].plot(x_d,y_d,'k',linewidth=lw)
				
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
				
				k+=1
	
	plt.savefig('../figures_resultats/chromosomes_collapsed'+str(ws)+'.svg') #coloring according to 3 quartet from original probabilities
	
	fig, axs = plt.subplots(ncols=5, nrows=4, sharey=True, sharex=True, figsize=(3.75, 10), layout="constrained")
	k = 1
	for row in range(4):
		for col in range(5):
			if k <=17:
				y_n = 1
				y_d = [0]
				x_d = [0]
				for i,j in enumerate(ret3[k][:-1]):
					y_d.append(j[0])
					x_d.append(-j[2]/d_max)
					y1 = j[0]
					y2 = (j[0]+ret3[k][i+1][0])/2
					y3 = ret3[k][i+1][0]
					c1 = cl[cl_3topos_dic[j[1]]] #collapse on 3 topos
					c2 = cl[cl_3topos_dic[ret3[k][i+1][1]]] #collapse on 3 topos
					if c1 != c2:
						axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [y2,y2],color=c1)
						y_n = y2
						
						axs[row,col].plot([0.0,1.0],[y2,y2],color='k',linewidth=lw)
				axs[row,col].fill_between([0.0,1.0],[ y_n,y_n], [ret3[k][-1][0],ret3[k][-1][0]],color=c1)
				axs[row,col].plot([0.0,1.0],[1,1],color='k',linewidth=lw)
				axs[row,col].plot([0.0,1.0],[ret3[k][-1][0],ret3[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([0.0,0.0],[1,ret3[k][-1][0]],color='k',linewidth=lw)
				axs[row,col].plot([1.0,1.0],[1,ret3[k][-1][0]],color='k',linewidth=lw)
				x_d.append(-ret3[k][-1][2]/d_max)
				y_d.append(ret3[k][-1][0])
				x_d.append(0)
				y_d.append(ret[k][-1][0])
				axs[row,col].plot(x_d,y_d,'k',linewidth=lw)
				
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
				
				k+=1
	
	plt.savefig('../figures_resultats/chromosomes_collapsed_sum'+str(ws)+'.svg') #coloring according to 3 quartet from the sum of probabilities

	
	fig, axs = plt.subplots(ncols=5, nrows=4, sharey=True, sharex=True, figsize=(7.5, 10), layout="constrained")
	k = 1
	for row in range(4):
		for col in range(5):
			if k < 18:
				minpos = 1
				maxpos = pval[k][-1][0]
				for i in range(1,16):
					axs[row,col].plot([i,i],[minpos,maxpos],color='blanchedalmond')
				for p in pval[k]:
					for i,j in enumerate(p[1]):
						pval_def = -np.exp(-lambda_pvals[i]*np.exp(j))+1
						b = binom(15,pval_def)
						pval_def_cor = 1-b.cdf(1)
# 						if j<=pvals10[i]:
						if pval_def_cor<=0.1:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='k')
							print(f"p-value on the window: {pval_def}")
							print(f"p-value on the window corrected for 15 tests: {pval_def_cor}")
				for p in pval[k]:
					for i,j in enumerate(p[1]):
						pval_def = -np.exp(-lambda_pvals[i]*np.exp(j))+1
						b = binom(15,pval_def)
						pval_def_cor = 1-b.cdf(1)
						if pval_def_cor<=0.05:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='r')
				for p in pval[k]:
					for i,j in enumerate(p[1]):
						pval_def = -np.exp(-lambda_pvals[i]*np.exp(j))+1
						b = binom(15,pval_def)
						pval_def_cor = 1-b.cdf(1)
						if pval_def_cor<=0.01:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='g')
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].tick_params(axis='y', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
			k+=1
	
	plt.savefig('../figures_resultats/probas_chromosomes_'+str(ws)+'.svg') # windows of likely gene transfers
