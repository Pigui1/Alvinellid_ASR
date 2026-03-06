import os
import numpy as np
from math import exp
import matplotlib.pyplot as plt
from random import sample
from copy import deepcopy
import pandas as pd
from subprocess import run, PIPE
from scipy.stats import binom
import multiprocessing as mp



def getASTRALscore():
	path = '/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/phylogenie/topologies/full_topos/'
	with open(path+'matrix_score.csv','w') as fout:
		fout.write('')
	with open(path+'matrix_quartet.csv','w') as fout:
		fout.write('')
	for i in range(1,16):
		tre1 = 'T'+str(i)+'.tre'
		sc1 = []
		qq1 = []
		for j in range(1,16):
			tre2 = 'T'+str(j)+'.tre'
		
			std = run('java -jar ~/Downloads/ASTRAL-master/Astral/astral.5.7.8.jar -i ' + path+tre1 + ' -q ' + path+tre2, stdout=PIPE, stderr=PIPE, shell=True, text=True)
			std = std.stderr.split('\n')
			for l in std:
				if 'Final quartet score is:' in l:
					score = l.split(': ')[1]
					sc1.append(score)
				if 'Number of quartet trees in the gene trees:' in l:
					quartets = l.split(': ')[1]
					qq1.append(quartets)
		with open(path+'matrix_score.csv','a') as fout:
			fout.write('\t'.join(sc1)+'\n')
		with open(path+'matrix_quartet.csv','a') as fout:
			fout.write('\t'.join(qq1)+'\n')



def get_duplicates_stats(r,ws,dis):
	ind = np.array(sample(range(len(r)),k=len(r)))
	
	pval_out = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] #p-values obtenues sur les fenêtres pour toutes les topologies
	
	for chromo in range(1,18):
		pos = 1
		max_pos = max(r['stop'].loc[chromo])
		
		mv = ws/10 #une fenêtre calculée tous les xxx nucléotides
		
		while pos <= max_pos:
			
			mask = (r['start'] >= pos-ws/2) & (r['stop'] <= pos+ws/2) & (r.index==chromo)
			ss = pd.DataFrame.sum(r.iloc[ind[mask]].iloc[:,22:37], axis=0)
			
# 			proba = np.matmul(ss,scores) #arbre coalescent
			proba = np.array(ss) #probabilités strictes
			
			b = binom(len(r.loc[mask]),dis)
			pvals = 1-b.cdf(np.array(proba))
			for i,p in enumerate(pvals):
				pval_out[i].append(p)
# 			print(chromo,pos,list(proba).index(max(proba))+1,st)
			
			pos += mv
			if pos < max_pos+mv:
				pos = min(pos,max_pos)
	
	return pval_out
			


def get_res(r,ws,dis):
	
	
	ret = {}
	pval = {}
	
	for chromo in range(1,18):
		ret[chromo] = []
		pval[chromo] = []
		phylo = -1
		pos = 1
		max_pos = max(r['stop'].loc[chromo])
		
		mv = ws/10 #une fenêtre calculée tous les xxx nucléotides
		
		while pos <= max_pos:
		
			mask = (r['start'] >= pos-ws/2) & (r['stop'] <= pos+ws/2) & (r.index==chromo)
			ss = pd.DataFrame.sum(r.loc[mask].iloc[:,22:37], axis=0)
			l_ali = pd.DataFrame.sum(r.loc[mask,'length'], axis=0)
# 			print(chromo,l_ali)
			
# 			proba = np.matmul(ss,scores) #arbre coalescent
			proba = np.array(ss) #probabilités strictes.
			
			ret[chromo].append([pos,list(proba).index(max(proba)),l_ali])
			
			b = binom(len(r.loc[mask]),dis[phylo])
			st = []
			pw = 1-b.cdf(ss)
# 			for i,p in enumerate(pw):
# 				if p <= max(pval[i],0.01): #0.01
# 					st.append(str(p))
# 				else:
# 					st.append('-')
			pval[chromo].append([pos,pw])
# 			st = '\t'.join(st)
# 			top = [str(chromo),str(pos),str(list(proba).index(max(proba))+1),st]
# 			print('\t'.join(top))
			
			pos += mv
			if pos < max_pos+mv:
				pos = min(pos,max_pos)
		
		
	return ret,pval
		
	
	
# 	print(r.iloc[ind[mask]]) #lot de gènes sue la fenêtre étudiée après rééchantillonnage 
# 	déplacer A
	


def get_prob(v):
	v2 = v-np.max(v)
	v2 = np.array([max(-30,x) for x in v2])
	v2 = np.exp(v2)
	return v2/np.sum(v2)



# getASTRALscore()



if __name__ == '__main__':
	p = 0.50 #filtre sur la probabilité minimale des gènes
	pr_al = 0.0 #filtre sur la proportion minimum de sites alignés par gène
	ws = 5e6 #taille de la fenêtre sur le génome
	dup0 = 1000 #nombre de réplicats aléatoires, minimum 20 pour p-value = 0.05, 100 pour p-value = 0.01


	r = []
	nbr_genes = 0
	with open('results_locus_probs.csv','r') as fin:
		l = fin.readline()
		col = l.strip('\n').split('\t')[1:]+['prop_ali']+['T'+str(i) for i in range(1,16)]
		ind = []
		for l in fin.readlines():
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
	

	# plt.hist(r['prop_ali'])
	# plt.show()
	# exit()


	# with open('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/phylogenie/topologies/full_topos/matrix_score.csv','r') as fin:
	# 	ll = fin.readlines()
	# scores = []
	# for l in ll:
	# 	l = l.strip('\n').split('\t')
	# 	scores.append([float(x) for x in l])
	# scores = np.matrix(scores)
	# 
	# with open('/Users/pgb/Documents/Data/articles/MDH_thermo_Alvinellidae/phylogenie/topologies/full_topos/matrix_quartet.csv','r') as fin:
	# 	ll = fin.readlines()
	# quartet = []
	# for l in ll:
	# 	l = l.strip('\n').split('\t')
	# 	quartet.append([float(x) for x in l])
	# quartet = np.matrix(quartet)
	# 
	# scores = scores/quartet
	# print(scores)


	r = r.loc[r['prop_ali']>pr_al]
	


	col = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']


	# p=0.5
	# q=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
	# qtot=[]
	# x=[]
	# while p<1.0:
	# 	qinter = []
	# 	for k in range(1,16):
	# 		L = 'L'+str(k)
	# 		s = r.loc[r[L]>p]
	# 		qinter.append(len(s))
	# 	for k in range(15):
	# 		q[k].append(qinter[k]/sum(qinter))
	# 	qtot.append(sum(qinter))
	# 	x.append(p)
	# 	p+=0.01
	# for i,qi in enumerate(q):
	# 	plt.plot(x,qi,color=col[i])
	# plt.plot(x,np.array(qtot),color='k')
	# plt.show()
	# exit()




	dis = []
	for k in range(1,16):
		L = 'L'+str(k)
		s = r.loc[r[L]>p]
		dis.append(len(s))
	dis = [x/sum(dis) for x in dis]


# 	ws = 1e6
# 	pvals_tot = []
# 	x_tot = []
# 
# 	while ws < 4e7:
# 		pval_tot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 
# 		dup = int(round(ws/1e6 * 500))
# 		with mp.Pool(processes=4) as pool:
# 			pvals_out = pool.starmap(get_duplicates_stats,[(r,ws,dis)]*dup)
# 		
# 		for pval in pvals_out:
# 			for i,_ in enumerate(pval_tot):
# 				pval_tot[i].append(min(pval[i])) #meilleure p-value obtenue par une fenêtre générée aléatoirement
# 		pvals = [0]*15
# 		for k,p in enumerate(pval_tot):
# 			p = sorted(p)
# 			pvals[k] = p[int(dup*0.0242-1)] #0.0242 : p-value 5% for 15 independant tested totpologies
# 		print(ws)
# 		x_tot.append(ws)
# 		pvals_tot.append(pvals)
# 		ws*=1.25
# 	pvals_tot = np.matrix(pvals_tot)
# 
# 	for i in range(15):
# 		plt.plot(x_tot,pvals_tot[:,i],color=col[i])
# 	plt.show()
# 
# 	exit()




	pval_tot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
	pval_tot_util = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# 	for k in range(dup):
	dup = int(round(ws/1e6 * dup0)) #adapter le nombre de simulations en fonction du nombre de fenêtres dans le génome
	with mp.Pool(processes=4) as pool:
		pvals_out = pool.starmap(get_duplicates_stats,[(r,ws,dis)]*dup)
	for pval in pvals_out:
		for i,_ in enumerate(pval_tot):
			pval_tot[i]+=pval[i]
			pval_tot_util[i].append(min(pval[i])) #meilleure p-value obtenue par une fenêtre générée aléatoirement
	
	
	pvals = [-1]*15
	for k,p in enumerate(pval_tot):
		p = sorted(p)
		x=[0]
		y=[0]
		for ind, i in enumerate(p[:-1]):
			if i != p[ind+1]:
				x.append(i)
				y.append((ind+1)/dup)
		plt.plot(x,y,color=col[k])
	plt.savefig('../figures_resultats/probas_per_window_'+str(ws)+'_2.svg')
# 	plt.show()
	
	for k,p in enumerate(pval_tot_util):
		p = sorted(p)
		pvals[k] = p[int(dup*0.0242-1)] #0.0242 : p-value 5% for 15 independant tested totpologies
	
	
	pvals1 = [-1]*15
	pvals10 = [-1]*15
	for k,p in enumerate(pval_tot_util):
		p = sorted(p)
		pvals1[k] = p[int(dup*0.01-1)] #0.01 : p-value 1% for 15 independant tested totpologies
		pvals10[k] = p[int(dup*0.036-1)] #0.036 : p-value 10% for 15 independant tested totpologies
		
	
	
	print(dup)
	print(pvals)
	
	with open('../figures_resultats/result_simulations_'+str(ws)+'.csv','w') as fout:
		for i in range(15):
			fout.write('\t'.join([str(x) for x in sorted(pval_tot_util[i])]) + '\n')
	
	
# 	pvals = [2.3e-5]*15 #5e6
# 	5e-6 #2.5e6
# 	5e-5 #4e7
	ret,pval = get_res(r,ws,dis)


	cl = ['#21618c', '#cb4335', '#117864', '#d35400', '#2471a3', '#239b56', '#e59866', '#138d75', '#3498db', '#5499c7', '#52be80', '#f39c12', '#48c9b0', '#a9cce3', '#f4d03f']

	
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
				x_d.append(0)
				y_d.append(ret[k][-1][0])
# 				axs[row,col].fill_between(x_d,y_d,color='k')#,[0,0],[1,ret[k][-1][0]])
				axs[row,col].plot(x_d,y_d,'k',linewidth=lw)
				
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
				
				k+=1
	
	plt.savefig('../figures_resultats/chromosomes_'+str(ws)+'_2.svg')
# 	plt.show()
	
	
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
						if j<=pvals10[i]:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='k')
				for p in pval[k]:
					for i,j in enumerate(p[1]):
						if j<=pvals[i]:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='r')
				for p in pval[k]:
					for i,j in enumerate(p[1]):
						if j<=pvals1[i]:
							axs[row,col].plot([i+1,i+1],[max(minpos,p[0]-ws/2),min(maxpos,p[0]+ws/2)],color='g')
				axs[row,col].tick_params(axis='x', which='both', bottom=False)
				axs[row,col].tick_params(axis='y', which='both', bottom=False)
				axs[row,col].set_frame_on(False)
			k+=1
	
	plt.savefig('../figures_resultats/probas_chromosomes_'+str(ws)+'_2.svg')
# 	plt.show()









