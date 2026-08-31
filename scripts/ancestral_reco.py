import ML_model
import matrices_model
from copy import deepcopy
import numpy as np
from scipy.linalg import expm



def calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,rightbranches,list_Anc,Anc_ter,mat,dic_eq,trif):
	#computing likelihood without root frequencies. call der=0 without derivative, der=1 with derivative for gamma, der=2 with derivative for matrix parameters
	#BL contains branch lengths and coordinates in the PCA of each mutation matrix
	
	for Anc in list_Anc:
		
		freq_1 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][0],eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
		freq_2 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][1],eigenvectors_1,eigenvectors_2,code)
		freq_3 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][2],eigenvectors_1,eigenvectors_2,code)
		freq_anc = freq_3
		
		
		mat_1 = {k:expm(BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(mat,freq_1,code)) for k in range(4)}
		mat_2 = {k:expm(BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(mat,freq_2,code)) for k in range(4)}
		if Anc == list_Anc[-1] and trif:
			mat_3 = {k:expm(BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(mat,freq_3,code)) for k in range(4)}
				
		for k in range(4):
			
			if rightbranches[Anc][0] == 1:
				probs_1 = np.matmul(mat_1[k],probs[k][tree[Anc][0]])
			else:
				probs_1 = np.matmul(probs[k][tree[Anc][0]].T,mat_1[k]).T
				freq_anc = freq_1
			if rightbranches[Anc][1] == 1:
				probs_2 = np.matmul(mat_2[k],probs[k][tree[Anc][1]])
			else:
				probs_2 = np.matmul(probs[k][tree[Anc][1]].T,mat_2[k]).T
				freq_anc = freq_2
			probs_tot = probs_1*probs_2
			
			if Anc == list_Anc[-1] and trif:
				if rightbranches[Anc][2] == 1:
					prob_3 = np.matmul(mat_3[k],probs[k][tree[Anc][2]])
				else:
					prob_3 = np.matmul(probs[k][tree[Anc][2]].T,mat_3[k]).T
				probs_tot *= prob_3			
			
			#in case the considered ancestor is the root of some positions
			if Anc == Anc_ter:
				compos = matrices_model.get_freq(compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code).reshape(-1,1)
			else:
				compos = freq_anc.reshape(-1,1)
			expo = np.log(compos) * roots[dic_eq[Anc]]
			prob_root = (np.e * np.ones(roots[dic_eq[Anc]].shape)) ** expo #multiplies by 1 if the residue is not a root, or by freq(residue) if it is the root
			probs_tot *= prob_root
			probs[k][Anc] = probs_tot
			
	return probs


def asr(Anc_cut,probs,tree,BL,descendants,Ancdes,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,mat):
	
	tree1, BL1, tree2, BL2, diceq = ML_model.cuttree(tree, BL, Anc_cut, descendants, Ancdes)
	Ancmax = 'Anc'+str(len(tree))
	
	if Anc_cut in Ancdes: #k n'est pas un descendant
		tree2_ancestors = {}
		for A in range(1,len(tree)+1):
			Anc = 'Anc'+str(A)
			tree2_ancestors[Anc] = [tree[Anc][0],tree[Anc][1]]
			if tree[Anc][0] in tree2_ancestors:
				tree2_ancestors[Anc] += tree2_ancestors[tree[Anc][0]]
			if tree[Anc][1] in tree2_ancestors:
				tree2_ancestors[Anc] += tree2_ancestors[tree[Anc][1]]
			if Anc == Anc_cut:
				break
		tree2_ancestors = set(tree2_ancestors[Anc_cut]+[Anc_cut])
		diceq2 = {diceq[key]:key for key in diceq if key not in tree2_ancestors} #equivalent ancêtres entre tree1 et tree
	else:
		diceq2 = {diceq[key]:key for key in diceq}

	rightbranch = {} #si l'orientation d'une branche est inversée entre treeb et t, alors dans le calcul des fréquences à l'équilibre on veut prendre les fréquences du fils sur cette branche

	for Anc in tree1:
		rightbranch[Anc] = [1,1,1]
		if tree1[Anc][0] in diceq2: #si le fils 0 est un ancêtre également
			if int(diceq2[Anc][3:]) < int(diceq2[tree1[Anc][0]][3:]): # Anc, dans l'ancien arbre, était le fils de tree1[Anc][0]
				rightbranch[Anc][0] = 0
		if tree1[Anc][1] in diceq2: #si le fils 1 est un ancêtre également
			if int(diceq2[Anc][3:]) < int(diceq2[tree1[Anc][1]][3:]): # Anc, dans l'ancien arbre, était le fils de tree1[Anc][1]
				rightbranch[Anc][1] = 0
		if tree1[Anc][2] in diceq2: #si le père est un ancêtre également
			if int(diceq2[Anc][3:]) > int(diceq2[tree1[Anc][2]][3:]): # Anc, dans l'ancien arbre, était le père de tree1[Anc][2]
				rightbranch[Anc][2] = 0


	Anc_ter = diceq['Anc'+str(len(tree))]
# 			anc_list = ['Anc'+str(x) for x in range(1,len(tree1)+1) if diceq2['Anc'+str(x)] not in node_done]
	anc_list = ['Anc'+str(x) for x in range(1,len(tree1)+1) if sum(rightbranch['Anc'+str(x)]) < 3]
	if len(anc_list) == 0: #anc_list is empy if we are at the root of the tree, and no branch needs to be reversed in the topology
		anc_list = ['Anc'+str(len(tree1))]
	probs_int = {k:{anc_loop:deepcopy(probs[k][diceq2[anc_loop]]) for anc_loop in tree1} for k in range(4)}
	for d in descendants:
		for k in range(4):
			probs_int[k][d] = deepcopy(probs[k][d])
	probs_int = calcul_probs_short(probs_int,tree1,BL1,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,rightbranch,anc_list,Anc_ter,mat,diceq2,0)
	p1 = {k:probs_int[k]['Anc'+str(len(tree1))] for k in range(4)}
	p2 = {k:probs[k][Anc_cut] for k in range(4)} #probability for second tree
	
	
	bl = BL1['Anc'+str(len(tree1))][2]
	freq = matrices_model.get_freq(compos_mean_1,compos_mean_2,bl,eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
	mat_w_f = matrices_model.equi_matrix(mat,freq,code)
	mat_equi = {k:expm(bl[0]*gam[k]*mat_w_f) for k in range(4)}
	
	prob_node = {}
	for k in range(4):
		prob_node[k] = p1[k]*np.matmul(mat_equi[k],p2[k])
	
	return prob_node


def calcul_probs(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param):
	#computing likelihood without root frequencies. call der=0 without derivative, der=1 with derivative for gamma, der=2 with derivative for matrix parameters
	#BL contains branch lengths and coordinates in the PCA of each mutation matrix
	
	if code == 'nt': #matrix without equilibrium frequencies
		mat = matrices_model.mat_gtr(*param,der)
	elif code == 'codon':
		mat = matrices_model.mat_codon(*param,der)
	elif code == 'aa':
		mat = model_aa
	elif code == 'codon_aa':
		mat = matrices_model.mat_codon(*param,der)
	
	
	for anc in range(1,len(tree)+1):
		Anc = 'Anc'+str(anc)
		
		freq_1 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][0],eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
		freq_2 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][1],eigenvectors_1,eigenvectors_2,code)
		freq_3 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][2],eigenvectors_1,eigenvectors_2,code)
		
		mat_1 = {k:expm(BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(mat,freq_1,code)) for k in range(4)}
		mat_2 = {k:expm(BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(mat,freq_2,code)) for k in range(4)}
		
		
		if anc == len(tree):
			mat_3 = {k:expm(BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(mat,freq_3,code)) for k in range(4)}
						
		for k in range(4):
			
			probs_1 = np.matmul(mat_1[k],probs[k][tree[Anc][0]])
			probs_2 = np.matmul(mat_2[k],probs[k][tree[Anc][1]])
			probs_tot = probs_1*probs_2
						
			if anc == len(tree):
				prob_3 = np.matmul(mat_3[k],probs[k][tree[Anc][2]])
				probs_tot *= prob_3			
			
			#in case the considered ancestor is the root of some positions
			if anc < len(tree):
				compos = freq_3.reshape(-1,1)
			else:
				compos = matrices_model.get_freq(compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code).reshape(-1,1)
			expo = np.log(compos) * roots[Anc]
			prob_root = (np.e * np.ones(roots[Anc].shape)) ** expo #multiplies by 1 if the residue is not a root, or by freq(residue) if it is the root
			probs_tot *= prob_root
			probs[k][Anc] = probs_tot
						
	return probs

