import numpy as np
from scipy.linalg import expm, expm_frechet
import PCA_module
import matrices_model
from scipy.integrate import quad
from scipy.special import gammainc
from math import gamma
from copy import deepcopy
from scipy.optimize import minimize
import multiprocessing as mp



def get_prob(seqs, tree, code): #builds matrices with probabilities of each state for the different species and nodes
	dic_out = {}
	
	codons = {'TTT':0,'TTC':1,'TTA':2,'TTG':3,'TCT':4,'TCC':5,'TCA':6,'TCG':7,
			  'TAT':8,'TAC':9,'TAA':10,'TAG':11,'TGT':12,'TGC':13,'TGA':14,'TGG':15,
			  'CTT':16,'CTC':17,'CTA':18,'CTG':19,'CCT':20,'CCC':21,'CCA':22,'CCG':23,
			  'CAT':24,'CAC':25,'CAA':26,'CAG':27,'CGT':28,'CGC':29,'CGA':30,'CGG':31,
			  'ATT':32,'ATC':33,'ATA':34,'ATG':35,'ACT':36,'ACC':37,'ACA':38,'ACG':39,
			  'AAT':40,'AAC':41,'AAA':42,'AAG':43,'AGT':44,'AGC':45,'AGA':46,'AGG':47,
			  'GTT':48,'GTC':49,'GTA':50,'GTG':51,'GCT':52,'GCC':53,'GCA':54,'GCG':55,
			  'GAT':56,'GAC':57,'GAA':58,'GAG':59,'GGT':60,'GGC':61,'GGA':62,'GGG':63}
	
	amino_acids = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
	
	nucleotides = {'A':0,'T':1,'C':2,'G':3,'A3':4,'T3':5,'C3':6,'G3':7}
	
	for s in seqs:
		if code == 'nt':
			dic_out[s] = np.zeros((8,len(seqs[s])))
			for i,c in enumerate(seqs[s]):
				if c == '-':
					if (i+1)%3 == 0:
						dic_out[s][:,i] = [0,0,0,0,1,1,1,1]
					else:
						dic_out[s][:,i] = [1,1,1,1,0,0,0,0]
				else:
					ind = ('3' if (i+1)%3 == 0 else '')
					dic_out[s][nucleotides[c+ind],i] = 1
		elif code == 'codon':
			max_i = int(len(seqs[s])/3)
			dic_out[s] = np.zeros((64,max_i))
			for i in range(max_i):
				if '-' in seqs[s][3*i:3*i+3]:
					dic_out[s][:,i] = [1]*64
				else:
					dic_out[s][codons[seqs[s][3*i:3*i+3]],i] = 1
		else:
			dic_out[s] = np.zeros((20,len(seqs[s])))
			for i,c in enumerate(seqs[s]):
				if c == '-':
					dic_out[s][:,i] = [1]*20
				else:
					dic_out[s][amino_acids[c],i] = 1
	
	for anc in tree:
		if code == 'nt':
			dic_out[anc] = np.zeros((8,len(seqs[s])))
		elif code == 'codon':
			dic_out[anc] = np.zeros((64,max_i))
		else:
			dic_out[anc] = np.zeros((20,len(seqs[s])))
	
	
	return dic_out


def get_prob_der(seqs, tree, code): #builds matrices with probabilities of each state for the different species and nodes
	dic_out = {}
	
	for s in seqs:
		if code == 'nt':
			dic_out[s] = np.zeros((8,len(seqs[s])))
		elif code == 'codon':
			dic_out[s] = np.zeros((64,max_i))
		else:
			dic_out[s] = np.zeros((20,len(seqs[s])))
	
	for anc in tree:
		if code == 'nt':
			dic_out[anc] = np.zeros((8,len(seqs[s])))
		elif code == 'codon':
			dic_out[anc] = np.zeros((64,max_i))
		else:
			dic_out[anc] = np.zeros((20,len(seqs[s])))
	
	return dic_out


def get_roots(seqs, tree, code):
	#identifying which ancestor is at the root, given gaps
	#rule: last ancestor which contains every contemporary descendant with residues
	#best to filter MSA to have all or almost all positions occupied by a residue, so that the root at a position is the true root of the tree (2 of the 3 branches from the root of the tree contain residues)
	
	Ancdes = {}
	for anc in range(1,len(tree)+1):
		Anc = 'Anc'+str(anc)
		Ancdes[Anc] = []
		if tree[Anc][0] not in tree:
			Ancdes[Anc].append(tree[Anc][0])
		else:
			Ancdes[Anc] += [x for x in Ancdes[tree[Anc][0]]]
		if tree[Anc][1] not in tree:
			Ancdes[Anc].append(tree[Anc][1])
		else:
			Ancdes[Anc] += [x for x in Ancdes[tree[Anc][1]]]
		if anc == len(tree):
			if tree[Anc][2] not in tree:
				Ancdes[Anc].append(tree[Anc][2])
			else:
				Ancdes[Anc] += [x for x in Ancdes[tree[Anc][2]]]
	
	ref_root = []
	if code == 'nt':
		root = {k:np.zeros((8,len(list(seqs.values())[0]))) for k in tree}
	if code == 'codon':
		root = {k:np.zeros((64,int(len(list(seqs.values())[0])/3))) for k in tree}
	else:
		root = {k:np.zeros((20,len(list(seqs.values())[0]))) for k in tree}
	
	if code == 'codon':
		for p in range(int(len(list(seqs.values())[0])/3)):
			nbr_res = 0
			for Des in Ancdes['Anc'+str(len(tree))]:
				if '-' not in seqs[Des][3*p:3*p+3]:
					nbr_res+=1
			for anc in range(1,len(tree)+1):
				Anc = 'Anc'+str(anc)
				k = 0
				for Des in Ancdes[Anc]:
					if '-' not in seqs[Des][3*p:3*p+3]:
						k+=1
				if nbr_res == k:
					root[Anc][:,p] += 1
					if anc < len(tree):
						ref_root.append((p,Anc))
					break
	else:
		for p in range(len(list(seqs.values())[0])):
			nbr_res = 0
			for Des in Ancdes['Anc'+str(len(tree))]:
				if '-' not in seqs[Des][p]:
					nbr_res+=1
			for anc in range(1,len(tree)+1):
				Anc = 'Anc'+str(anc)
				k = 0
				for Des in Ancdes[Anc]:
					if '-' != seqs[Des][p]:
						k+=1
				if nbr_res == k:
					root[Anc][:,p] += 1
					if anc<len(tree):
						ref_root.append((p,Anc))
					break
	
	return root,ref_root


def gammainc2(x, alpha, beta):
	return gammainc(alpha, beta*x)


def gammaden(x, alpha, beta):
	return ((x*beta)**alpha) * np.exp(-beta*x) / gamma(alpha)


def classesgamma(g): #obtenir la moyenne du taux d'évolution selon la loi gamma pour 5 classes (0, 0.2, 0.4, 0.6, 0.8, 1)
	nbrclasses = 4
	classes = [i*1/nbrclasses for i in range(1,nbrclasses)]
	bornesout = []
	alpha = g
	beta = g
	for k in classes:
		binf = -10
		bsup = 10 #arbitraire, vérifié en dessous
		while gammainc2(bsup,alpha,beta) < k: #fonction de distribution de la loi gamma
			bsup = bsup*2
		while gammainc2(binf,alpha,beta) > k:
			binf = binf*2
		t = (binf+bsup)/2
		while abs(gammainc2(t,alpha,beta)-k) > 10e-15: #seuil arbitraire
			if gammainc2(t,alpha,beta) > k:
				bsup = t
			elif gammainc2(t,alpha,beta) < k:
				binf = t
			t = (binf+bsup)/2
		bornesout.append(t)
	
	classesout = [0]*(len(classes)+1)
	
	for b, _ in enumerate(bornesout[:-1]):
		classesout[b+1] = (len(classes)+1)*quad(gammaden, bornesout[b], bornesout[b+1], args=(alpha, beta))[0]
	
	classesout[0] = (len(classes)+1)*quad(gammaden, 0, bornesout[0], args=(alpha, beta))[0]
	classesout[-1] = (len(classes)+1)*quad(gammaden, bornesout[-1], np.inf, args=(alpha, beta))[0]
	
	return np.array(classesout)


def classesgamma_der(g):
	eps = 1e-5
	vec_1 = classesgamma(g-eps)
	vec_2 = classesgamma(g+eps)
	return (vec_2-vec_1)/(2*eps)


def calcul_probs(probs,tree,BL,roots,coord_root,g_i,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param,der):
	#computing likelihood without root frequencies. call der=0 without derivative, der=1 with derivative for gamma, der=2 with derivative for matrix parameters
	#BL contains branch lengths and coordinates in the PCA of each mutation matrix
	
	gam = classesgamma(g_i)
	
	if der < 2:
		if code == 'nt': #matrix without equilibrium frequencies
			mat = matrices_model.mat_gtr(*param,der)
		elif code == 'codon':
			mat = matrices_model.mat_codon(*param,der)
		elif code == 'aa':
			mat = model_aa
		elif code == 'codon_aa':
			mat = matrices_model.mat_codon(*param,der)
	
	elif der == 2:
		if code == 'nt': #matrix without equilibrium frequencies
			mat,mat_der = matrices_model.mat_gtr(*param,der)
		elif code == 'codon':
			mat,mat_der = matrices_model.mat_codon(*param,der)
		elif code == 'codon_aa':
			mat,mat_der = matrices_model.mat_codon(*param,der)
	
	probs_der = {}
	for k in range(4):
		if der == 1:
			probs_der[k] = {x:np.zeros(probs[k][x].shape) for x in probs[k]}
		if der == 2:
			probs_der[k] = {x:np.zeros((len(param),*probs[k][x].shape)) for x in probs[k]}
	
	
	for anc in range(1,len(tree)+1):
		Anc = 'Anc'+str(anc)
		
		freq_1 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][0],eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
		freq_2 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][1],eigenvectors_1,eigenvectors_2,code)
		freq_3 = matrices_model.get_freq(compos_mean_1,compos_mean_2,BL[Anc][2],eigenvectors_1,eigenvectors_2,code)
		
		mat_1 = {k:expm(BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(mat,freq_1,code)) for k in range(4)}
		mat_2 = {k:expm(BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(mat,freq_2,code)) for k in range(4)}
		
		
		if anc == len(tree):
			mat_3 = {k:expm(BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(mat,freq_3,code)) for k in range(4)}
		if der == 1:
			g_der = classesgamma_der(g_i)
			mat_1_der = {k:(BL[Anc][0][0]*g_der[k]*matrices_model.equi_matrix(mat,freq_1,code)).dot(expm(BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(mat,freq_1,code))) for k in range(4)} #derivative according to gamma
			mat_2_der = {k:(BL[Anc][1][0]*g_der[k]*matrices_model.equi_matrix(mat,freq_2,code)).dot(expm(BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(mat,freq_2,code))) for k in range(4)} #derivative according to gamma
			if anc == len(tree):
				mat_3_der = {k:(BL[Anc][2][0]*g_der[k]*matrices_model.equi_matrix(mat,freq_3,code)).dot(expm(BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(mat,freq_3,code))) for k in range(4)} #derivative according to gamma
		if der == 2:
			#ATTENTION erreur avec equilibrage des frequences * dd, ca ne fonctionne pas comme ça. peut être que la fonction equi_matrix fonctionnerait
			#Surement erreur avec la normalisation des matrices, voir fonction L_branch
			mat_1_der = {k:[expm_frechet(BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(mat,freq_1,code), BL[Anc][0][0]*gam[k]*matrices_model.equi_matrix(dd,freq_1,code), compute_expm=False, check_finite=False) for dd in mat_der] for k in range(4)} #derivative according to matrix parameters
			mat_2_der = {k:[expm_frechet(BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(mat,freq_2,code), BL[Anc][1][0]*gam[k]*matrices_model.equi_matrix(dd,freq_2,code), compute_expm=False, check_finite=False) for dd in mat_der] for k in range(4)} #derivative according to matrix parameters
			if anc == len(tree):
				mat_3_der = {k:[expm_frechet(BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(mat,freq_3,code), BL[Anc][2][0]*gam[k]*matrices_model.equi_matrix(dd,freq_3,code), compute_expm=False, check_finite=False) for dd in mat_der] for k in range(4)} #derivative according to matrix parameters
		
		for k in range(4):
			
			probs_1 = np.matmul(mat_1[k],probs[k][tree[Anc][0]])
			probs_2 = np.matmul(mat_2[k],probs[k][tree[Anc][1]])
			probs_tot = probs_1*probs_2
			
			if der==1:
				probs_1_der = np.matmul(mat_1_der[k],probs[k][tree[Anc][0]]) + np.matmul(mat_1[k],probs_der[k][tree[Anc][0]])
				probs_2_der = np.matmul(mat_2_der[k],probs[k][tree[Anc][1]]) + np.matmul(mat_2[k],probs_der[k][tree[Anc][1]])
				probs_tot_der = probs_1_der*probs_2 + probs_1*probs_2_der
				
			if der==2:
				probs_1_der = [np.matmul(mat_1_der[k][x],probs[k][tree[Anc][0]]) + np.matmul(mat_1[k],probs_der[k][tree[Anc][0]][x]) for x in range(len(param))]
				probs_2_der = [np.matmul(mat_2_der[k][x],probs[k][tree[Anc][1]]) + np.matmul(mat_2[k],probs_der[k][tree[Anc][1]][x]) for x in range(len(param))]
				probs_tot_der = [probs_1_der[x]*probs_2 + probs_1*probs_2_der[x] for x in range(len(param))]
			
			if anc == len(tree):
				prob_3 = np.matmul(mat_3[k],probs[k][tree[Anc][2]])
				probs_tot *= prob_3
				if der==1:
					prob_3_der = np.matmul(mat_3_der[k],probs[k][tree[Anc][2]]) + np.matmul(mat_3[k],probs_der[k][tree[Anc][2]])
					probs_tot_der = probs_tot_der*prob_3 + probs_1*probs_2*prob_3_der
				if der==2:
					probs_3_der = [np.matmul(mat_3_der[k][x],probs[k][tree[Anc][2]]) + np.matmul(mat_3[k],probs_der[k][tree[Anc][2]][x]) for x in range(len(param))]
					probs_tot_der = [probs_tot_der[x]*prob_3 + probs_1*probs_2*prob_3_der[x] for x in range(len(param))]
			
			
			#in case the considered ancestor is the root of some positions
			if anc < len(tree):
				compos = freq_3.reshape(-1,1)
			else:
				compos = matrices_model.get_freq(compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code).reshape(-1,1)
			expo = np.log(compos) * roots[Anc]
			prob_root = (np.e * np.ones(roots[Anc].shape)) ** expo #multiplies by 1 if the residue is not a root, or by freq(residue) if it is the root
			probs_tot *= prob_root
			probs[k][Anc] = probs_tot
			
			if der == 1:
				probs_tot_der *= prob_root
				probs_der[k][Anc] = probs_tot_der
			if der == 2:
				probs_tot_der = [x*prob_root for x in probs_tot_der]
				probs_der[k][Anc] = probs_tot_der
			
	
	if der==1 or der==2:
		return probs, probs_der
	else:
		return probs


def calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranches,list_Anc,Anc_ter,mat,dic_eq,trif):
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


def iter_branch_parall(Anc_cut,probs,tree,BL,descendants,Ancdes,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,mat):
	
	tree1, BL1, tree2, BL2, diceq = cuttree(tree, BL, Anc_cut, descendants, Ancdes)
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
	probs_int = calcul_probs_short(probs_int,tree1,BL1,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranch,anc_list,Anc_ter,mat,diceq2,0)
	p1 = {k:probs_int[k]['Anc'+str(len(tree1))] for k in range(4)}
	p2 = {k:probs[k][Anc_cut] for k in range(4)} #probability for second tree
	
	bl_res = minimize(L_branch,BL1['Anc'+str(len(tree1))][2],(p1,p2,mat,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code),method='TNC',jac=True,bounds=[(1e-6,np.inf)]+[(-np.inf,np.inf) for x in BL1['Anc'+str(len(tree1))][2][1:]])
	bl_res = bl_res.x
	
	L_branch_1 = L_branch_simple(BL1['Anc'+str(len(tree1))][2],p1,p2,mat,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code)
	L_branch_2 = L_branch_simple(bl_res,p1,p2,mat,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code)
	
	if (L_branch_1 - L_branch_2)**2 < 1e-4: #if the L after that branch optimization didn't change, then the probability around Anc_cut won't change further and doesn't need to be calculated
		return [Anc_cut,bl_res]
	else:
		return [0,bl_res]


def iter_branch(L_init,probs,tree,BL,descendants,roots,coord_root,g_i,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,param,code,model_aa,ref_root,nt):
	node_iter = descendants + ['Anc'+str(i) for i in range(1,len(tree))] #nodes that need to be optimized
	node_done = []
	
	#calcul of the fixed parameters
	gam = classesgamma(g_i)
	
	if code == 'nt': #matrix without equilibrium frequencies
		mat = matrices_model.mat_gtr(*param,0)
	elif code == 'codon':
		mat = matrices_model.mat_codon(*param,0)
	elif code == 'aa':
		mat = model_aa
	elif code == 'codon_aa':
		mat = matrices_model.mat_codon(*param,0)
	
	#calcul of the initial likelihood, before branch and coefficients optimization
	Ancmax = 'Anc'+str(len(tree))
	
	Ancdes = {}
	for anc in range(1,len(tree)+1):
		Anc = 'Anc'+str(anc)
		Ancdes[Anc] = []
		if tree[Anc][0] not in tree:
			Ancdes[Anc].append(tree[Anc][0])
		else:
			Ancdes[Anc] += [x for x in Ancdes[tree[Anc][0]]]
		if tree[Anc][1] not in tree:
			Ancdes[Anc].append(tree[Anc][1])
		else:
			Ancdes[Anc] += [x for x in Ancdes[tree[Anc][1]]]
		if anc == len(tree):
			if tree[Anc][2] not in tree:
				Ancdes[Anc].append(tree[Anc][2])
			else:
				Ancdes[Anc] += [x for x in Ancdes[tree[Anc][2]]]
	
	
	L = 0
	nt = min(nt,4) #we are not calculating more than 4 branches at the same time as the performance decreases
	while (L-L_init)**2 > 1e-4:
		L = L_init
		node_iter = [x for x in node_iter if x not in node_done]
		i = 0
		param1 = []
		while i < len(node_iter):
			param1.append(node_iter[i:i+nt])
			i+=nt
		
		for k in param1:
			param = [(x,probs,tree,BL,descendants,Ancdes,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,mat) for x in k]
		
			with mp.Pool(processes = nt) as pool:
				res = pool.starmap(iter_branch_parall,param)
			pool.close()
			
			node_done1 = [x[0] for x in res if x != 0]
			node_done += node_done1
			bl_res_tot = [x[1] for x in res]
			
			for Anc_cut, bl_res in zip(k,bl_res_tot):
				for anc_it in tree:
					if Anc_cut in tree[anc_it][0:2]:
						BL[anc_it][tree[anc_it].index(Anc_cut)] = [x for x in bl_res]
					if anc_it == Ancmax and Anc_cut == tree[anc_it][2]:
						BL[anc_it][2] = [x for x in bl_res]
				if Anc_cut in tree:
					BL[Anc_cut][2] = [x for x in bl_res]
			
			
			anc_to_calc_tot = []
			rightbranches = {k:[1,1,1] for k in BL}
			dic_eq = {k:k for k in tree}
			
			for Anc_cut in k:
				Anc2 = Anc_cut
				if Anc_cut not in tree:
					for t in tree:
						if Anc_cut in tree[t]:
							anc_to_calc = [t]
							Anc2 = t
				else:
					anc_to_calc = []
				while Anc2 != Ancmax:
					anc_to_calc.append(tree[Anc2][2])
					Anc2 = tree[Anc2][2]
				anc_to_calc_tot += anc_to_calc
			anc_to_calc_tot = set(anc_to_calc)
			anc_to_calc = ['Anc'+str(x) for x in range(1,len(tree)+1) if 'Anc'+str(x) in anc_to_calc_tot] #only updating nodes that are changed by the new branch lengths
			
			probs = calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranches,anc_to_calc,Ancmax,mat,dic_eq,1)
		
		if len(coord_root)>1:
# 			diceq_n = {x:x for x in tree}
# 			rightbranch = {x:[1,1,1] for x in BL}
# 			probs = calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranch,['Anc'+str(x) for x in range(1,len(tree)+1)],Ancmax,mat,diceq_n,1)
			compos = matrices_model.get_freq(compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code).reshape(-1,1)
			expo = np.log(compos) * roots[Ancmax]
			prob_root = (np.e * np.ones(roots[Ancmax].shape)) ** expo #multiplies by 1 if the residue is not a root, or by freq(residue) if it is the root
			p1 = {k:probs[k][Ancmax]/prob_root for k in range(4)}
	
		# 	L, L_der = optim_root(roots[Ancmax],p1,ref_root,compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code)
			root_res = minimize(optim_root,coord_root,(roots[Ancmax],p1,ref_root,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code),method='TNC',jac=True)
			root_res = root_res.x
			coord_root = [x for x in root_res]
		
		
		rightbranches = {k:[1,1,1] for k in BL}
		dic_eq = {k:k for k in tree}
# 		probs = calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranches,['Anc'+str(k) for k in range(1,len(tree)+1)],Ancmax,mat,dic_eq,1)
		probs = calcul_probs_short(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,rightbranches,[Ancmax],Ancmax,mat,dic_eq,1)
		L_init = np.zeros((1,probs[0][Ancmax].shape[1]))
		for k in range(4):
			for i,a in ref_root:
				probs[k][Ancmax][:,i] = probs[k][a][:,i]
			L_init+= np.sum(probs[k][Ancmax],axis=0)
		L_init = -np.sum(np.log(L_init/4))
		
		print(L_init)
	
	return BL, coord_root
	

def optim_root(coord_root,roots,p1,ref_root,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code):
	
	compos = matrices_model.get_freq(compos_mean_1,compos_mean_2,coord_root,eigenvectors_1,eigenvectors_2,code).reshape(-1,1)
	expo = np.log(compos) * roots
	prob_root = (np.e * np.ones(roots.shape)) ** expo #multiplies by 1 if the residue is not a root, or by freq(residue) if it is the root
	p = {k:p1[k]*prob_root for k in range(4)}
	
	if min(compos) > 0.00015:
		compos_der = eigenvectors_1
	else:
		eps = 1e-10
		coord_pert = [[x for x in coord_root]for i in range(len(eigenvectors_1))]
		for i in range(len(eigenvectors_1)):
			coord_pert[i][i+1] = coord_pert[i][i+1] + eps
		compos_2 = [matrices_model.get_freq(compos_mean_1,compos_mean_2,crp,eigenvectors_1,eigenvectors_2,code) for crp in coord_pert]
		compos_der = [(c2-compos.reshape(1,-1))/eps for c2 in compos_2]
	
	expo_der = [(x.reshape(-1,1)/compos) * roots for x in compos_der]
	prob_root_der = [x * prob_root for x in expo_der]
	p_der = [{k:p1[k]*x for k in range(4)} for x in prob_root_der]	
	
	L = np.zeros((1,p[0].shape[1]))
	L_der = np.zeros((len(p_der),1,p[0].shape[1]))
	for k in range(4):
		for i,a in ref_root:
			p[k][:,i] = np.ones(p[0].shape[0])
			for der in range(len(p_der)):
				p_der[der][k][:,i] = np.zeros(p[0].shape[0])
		Lk = np.sum(p[k],axis=0)
		L += Lk
		for der in range(len(p_der)):
			L_der[der] += np.sum(p_der[der][k],axis=0)
	
	L_der = [-np.sum(der/L) for der in L_der]
	L = -np.sum(np.log(L/4))
	
	return L, [0]+L_der


def L_branch(bl,p1,p2,mat,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code):
	freq = matrices_model.get_freq(compos_mean_1,compos_mean_2,bl,eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
	mat_w_f = matrices_model.equi_matrix(mat,freq,code)
	mat_equi = {k:expm(bl[0]*gam[k]*mat_w_f) for k in range(4)}
	
	mat_equi_der_t = {k:np.matmul(gam[k]*mat_w_f,mat_equi[k]) for k in range(4)}
	
	freq_der = eigenvectors_1
	
	if min(freq) > 0.00015: #freq has been rescaled to avoid negative frequencies
		freq_der = eigenvectors_1
	else:
		eps = 1e-10
		bl_pert = [[x for x in bl]for i in range(len(eigenvectors_1))]
		for i in range(len(eigenvectors_1)):
			bl_pert[i][i+1] = bl_pert[i][i+1] + eps	
		freq_2 = [matrices_model.get_freq(compos_mean_1,compos_mean_2,blp,eigenvectors_1,eigenvectors_2,code) for blp in bl_pert]
		freq_der = [(f2-freq)/eps for f2 in freq_2]
	
	
	alp = np.sum(freq*((freq*mat).T))
	alp_der = [np.sum(fd*((freq*mat).T) + freq*((fd*mat).T)) for fd in freq_der]
	
	
	directions = [(fd/alp - freq*alpp/(alp**2))*mat for fd,alpp in zip(freq_der,alp_der)]
	for i in range(len(directions)):
		directions[i] = directions[i] - np.sum(directions[i],axis=1)*np.identity(20)
	
	mat_der_coeff = [{k:expm_frechet(bl[0]*gam[k]*mat_w_f, bl[0]*gam[k]*dr, compute_expm=False, check_finite=False) for k in range(4)} for dr in directions]
	
	L = np.zeros((1,p1[0].shape[1]))
	L_der_t = np.zeros((1,p1[0].shape[1]))
	L_der_coeff = np.zeros((len(mat_der_coeff),1,p1[0].shape[1]))
	for k in range(4):
		probs = p1[k]*np.matmul(mat_equi[k],p2[k])
		probs_der_t = p1[k]*np.matmul(mat_equi_der_t[k],p2[k])
		Lk = np.sum(probs,axis=0)
		L_derk = np.sum(probs_der_t,axis=0)
		L += Lk
		L_der_t += L_derk
		for der in range(len(mat_der_coeff)):
			probs_der = p1[k]*np.matmul(mat_der_coeff[der][k],p2[k])
			L_der_coeff[der] += np.sum(probs_der,axis=0)
	
	L_der_t = -np.sum(L_der_t/L)
	L_der_coeff = [-np.sum(der/L) for der in L_der_coeff]
	L = -np.sum(np.log(L/4))
	
	return L, [L_der_t] + L_der_coeff


def L_branch_simple(bl,p1,p2,mat,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code):
	freq = matrices_model.get_freq(compos_mean_1,compos_mean_2,bl,eigenvectors_1,eigenvectors_2,code) #equilibrium frequencies around the node
	mat_w_f = matrices_model.equi_matrix(mat,freq,code)
	mat_equi = {k:expm(bl[0]*gam[k]*mat_w_f) for k in range(4)}
	
	L = np.zeros((1,p1[0].shape[1]))
	for k in range(4):
		probs = p1[k]*np.matmul(mat_equi[k],p2[k])
		Lk = np.sum(probs,axis=0)
		L += Lk
	L = -np.sum(np.log(L/4))
	
	return L	


def cuttree(tree, BL, Anceval, Desc, Ancdes): #reracine un arbre selon le noeud parent de Anceval
	
	tmaxinit = len(tree)
	if not (Anceval == 'Anc'+str(tmaxinit) or Anceval == tree['Anc'+str(tmaxinit)][2]):
		
		if Anceval not in Desc:
			iAnceval = int(Anceval[3:])
			tmax = int(len(tree) - len(Ancdes[Anceval]) +1)
			diceq = {tree[Anceval][2]:'Anc'+str(tmax)}	
			inverse = [Anceval]
			Ancevalparent = tree[Anceval][2]
		else:
			tmax = int(len(tree))
			inverse = [Anceval]
			for Anc in tree:
				if Anceval in tree[Anc]:
					iAnceval = int(Anc[3:])
					Ancevalparent = Anc
					break
			diceq = {'Anc'+str(iAnceval):'Anc'+str(tmax)}	
			inverse = [Anceval, 'Anc'+str(iAnceval)]
		while iAnceval != tmaxinit and 'Anc'+str(iAnceval) != tree['Anc'+str(tmaxinit)][2]: #peut-être un soucis ici, pas testé dans tous les cas
			inverse.append(tree['Anc'+str(iAnceval)][2])
			iAnceval = int(tree['Anc'+str(iAnceval)][2][3:])
	
	
		i = tmax
		explore = []
		if tree[Ancevalparent][0] == Anceval: #si le fils 0 est Anceval, c'est cette branche qui passe en ancestrale
			if tree[Ancevalparent][1] not in Desc: #si le fils 1 est un noeud interne
				i-=1
				son0 = 'Anc'+str(i) #on le nomme 'Anc' décrémenté
				explore.append(tree[Ancevalparent][1])
				diceq[tree[Ancevalparent][1]] = 'Anc'+str(i)
			else:
				son0 = tree[Ancevalparent][1] #on lui donne le nom de l'espèce
			BL0 = BL[Ancevalparent][1]
			BL2 = BL[Ancevalparent][0]
		elif tree[Ancevalparent][1] == Anceval: #si le fils 1 est Anceval, c'est cette branche qui passe en ancestrale
			if tree[Ancevalparent][0] not in Desc: #si le fils 1 est un noeud interne
				i-=1
				son0 = 'Anc'+str(i) #on le nomme 'Anc' décrémenté
				explore.append(tree[Ancevalparent][0])
				diceq[tree[Ancevalparent][0]] = 'Anc'+str(i)
			else:
				son0 = tree[Ancevalparent][0] #on lui donne le nom de l'espèce
			BL0 = BL[Ancevalparent][0]
			BL2 = BL[Ancevalparent][1]
		if tree[Ancevalparent][2] not in Desc: #si l'ancêtre est un noeud interne
			i-=1
			son1 = 'Anc'+str(i) #on le nomme 'Anc' décrémenté
			explore.append(tree[Ancevalparent][2])
			diceq[tree[Ancevalparent][2]] = 'Anc'+str(i)
		else:
			son1 = tree[Ancevalparent][2] #on lui donne le nom de l'espèce
		BL1 = BL[Ancevalparent][2]
		treeout = {'Anc'+str(tmax):[son0,son1,'XYZ']}
		BLout = {'Anc'+str(tmax):[BL0,BL1,BL2]}
		
		while len(explore)>0: #i>1
			explore2 = []
			for Ancexp in explore:
				if Ancexp in inverse: #on est sur le chemin entre l'ancêtre absolu de l'arbre et l'ancêtre qu'on veut mettre à la racine, il va falloir faire pivoter la topologie
					if tree[Ancexp][0] not in diceq: #si le fils 0 de l'ancêtre qu'on explore n'est pas le nouvel ancêtre
						if tree[Ancexp][0] not in Desc:
							i-=1
							son0 = 'Anc'+str(i)
							explore2.append(tree[Ancexp][0])
							diceq[tree[Ancexp][0]] = 'Anc'+str(i)
							BL0 = BL[Ancexp][0]
						else:
							son0 = tree[Ancexp][0]
							BL0 = BL[Ancexp][0]
						son2 = diceq[tree[Ancexp][1]]
						BL2 = BL[Ancexp][1]
					elif tree[Ancexp][1] not in diceq: #si le fils 1 de l'ancêtre qu'on explore n'est pas le nouvel ancêtre
						if tree[Ancexp][1] not in Desc:
							i-=1
							son0 = 'Anc'+str(i)
							explore2.append(tree[Ancexp][1])
							diceq[tree[Ancexp][1]] = 'Anc'+str(i)
							BL0 = BL[Ancexp][1]
						else:
							son0 = tree[Ancexp][1]
							BL0 = BL[Ancexp][1]
						son2 = diceq[tree[Ancexp][0]]
						BL2 = BL[Ancexp][0]
					if tree[Ancexp][2] not in Desc: #l'ancêtre de tree devient le deuxième fils de treeout
						i-=1
						son1 = 'Anc'+str(i)
						explore2.append(tree[Ancexp][2])
						diceq[tree[Ancexp][2]] = 'Anc'+str(i)
						BL1 = BL[Ancexp][2]
					else:
						son1 = tree[Ancexp][2]	
						BL1 = BL[Ancexp][2]			
				
				else:
					if tree[Ancexp][0] not in Desc:
						i-=1
						son0 = 'Anc'+str(i)
						explore2.append(tree[Ancexp][0])
						diceq[tree[Ancexp][0]] = 'Anc'+str(i)
						BL0 = BL[Ancexp][0]
					else:
						son0 = tree[Ancexp][0]
						BL0 = BL[Ancexp][0]
					if tree[Ancexp][1] not in Desc:
						i-=1
						son1 = 'Anc'+str(i)
						explore2.append(tree[Ancexp][1])
						diceq[tree[Ancexp][1]] = 'Anc'+str(i)
						BL1 = BL[Ancexp][1]
					else:
						son1 = tree[Ancexp][1]
						BL1 = BL[Ancexp][1]
					son2 = diceq[tree[Ancexp][2]]
					BL2 = BL[Ancexp][2]
				
				treeout[diceq[Ancexp]] = [son0, son1, son2]
				BLout[diceq[Ancexp]] = [BL0, BL1, BL2]
			explore = [j for j in explore2]
	
	elif Anceval == 'Anc'+str(tmaxinit) or Anceval == tree['Anc'+str(tmaxinit)][2]:
		
		Ancevalentry = Anceval
		if Anceval in Desc:
			Anceval = 'Anc'+str(tmaxinit)
		else:
			Anceval = tree[Anceval][2]
		if Anceval not in Desc:
			tmax = len(Ancdes[Anceval]) - 2 
			i = tmax
			explore = []
			diceq = {Anceval : 'Anc'+str(tmax)}
	
			if tree[Anceval][0] in Desc:
				son0 = tree[Anceval][0]
			else:
				i-=1
				son0 = 'Anc' +str(i)
				explore.append(tree[Anceval][0])
				diceq[tree[Anceval][0]] = son0
			if tree[Anceval][1] in Desc:
				son1 = tree[Anceval][1]
			else:
				i-=1
				son1 = 'Anc' +str(i)
				explore.append(tree[Anceval][1])
				diceq[tree[Anceval][1]] = son1
			BL0 = BL[Anceval][0]
			BL1 = BL[Anceval][1]
			BL2 = BL[Anceval][2]
			treeout = {'Anc'+str(tmax) : [son0,son1,'XYZ']}
			BLout = {'Anc'+str(tmax) : [BL0,BL1,BL2]}
		
			while len(explore)>0:
				explore2 = []
				for Ancexp in explore:
					if tree[Ancexp][0] in Desc:
						son0 = tree[Ancexp][0]
					else:
						i-=1
						son0 = 'Anc' +str(i)
						explore2.append(tree[Ancexp][0])
						diceq[tree[Ancexp][0]] = son0
					if tree[Ancexp][1] in Desc:
						son1 = tree[Ancexp][1]
					else:
						i-=1
						son1 = 'Anc' +str(i)
						explore2.append(tree[Ancexp][1])
						diceq[tree[Ancexp][1]] = son1
					BL0 = BL[Ancexp][0]
					BL1 = BL[Ancexp][1]
					BL2 = BL[Ancexp][2]
					son2 = diceq[tree[Ancexp][2]]
					treeout[diceq[Ancexp]] = [son0,son1,son2]
					BLout[diceq[Ancexp]] = [BL0,BL1,BL2]
				explore = [j for j in explore2]
		else:
			treeout = {}
			BLout = {}
		Anceval = Ancevalentry

	
	#headache
		
	if Anceval not in Desc:
		tmax = len(Ancdes[Anceval]) - 1
		i = tmax
		explore = []
		diceq[Anceval] = 'Anc'+str(tmax)
	
		if tree[Anceval][0] in Desc:
			son0 = tree[Anceval][0]
		else:
			i-=1
			son0 = 'Anc' +str(i)
			explore.append(tree[Anceval][0])
			diceq[tree[Anceval][0]] = son0
		if tree[Anceval][1] in Desc:
			son1 = tree[Anceval][1]
		else:
			i-=1
			son1 = 'Anc' +str(i)
			explore.append(tree[Anceval][1])
			diceq[tree[Anceval][1]] = son1
		BL0 = BL[Anceval][0]
		BL1 = BL[Anceval][1]
		BL2 = BL[Anceval][2]
		treeout2 = {'Anc'+str(tmax) : [son0,son1,'XYZ']}
		BLout2 = {'Anc'+str(tmax) : [BL0,BL1,BL2]}
		
		while len(explore)>0:
			explore2 = []
			for Ancexp in explore:
				if tree[Ancexp][0] in Desc:
					son0 = tree[Ancexp][0]
				else:
					i-=1
					son0 = 'Anc' +str(i)
					explore2.append(tree[Ancexp][0])
					diceq[tree[Ancexp][0]] = son0
				if tree[Ancexp][1] in Desc:
					son1 = tree[Ancexp][1]
				else:
					i-=1
					son1 = 'Anc' +str(i)
					explore2.append(tree[Ancexp][1])
					diceq[tree[Ancexp][1]] = son1
				BL0 = BL[Ancexp][0]
				BL1 = BL[Ancexp][1]
				BL2 = BL[Ancexp][2]
				son2 = diceq[tree[Ancexp][2]]
				treeout2[diceq[Ancexp]] = [son0,son1,son2]
				BLout2[diceq[Ancexp]] = [BL0,BL1,BL2]
			explore = [j for j in explore2]
	else:
		treeout2 = {}
		BLout2 = {}
	#deuxième arbre:
	#prendre Anceval
	#compter combien de noeuds avec Ancdes
	#renommer les Ancetres en descendant en partant de Anceval
	#on recopie exactement l'arbre en changeant simplement les noms avec un dictionnaire d'équivalences
	#BL est identique	
	
	return treeout, BLout, treeout2, BLout2, diceq


def calcul_L(gam,probs,tree,BL,roots,coord_root,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param,ref_root,der):
	Ancmax = 'Anc'+str(len(tree))
	if der == 0:
		probs = calcul_probs(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param,der)
# 		print(probs[0])
		L = np.zeros((1,probs[0][Ancmax].shape[1]))
		for k in range(4):
			for i,a in ref_root:
				probs[k][Ancmax][:,i] = probs[k][a][:,i]
			L+= np.sum(probs[k][Ancmax],axis=0)
		L = np.sum(np.log(L/4))
		
		return -L
	
	if der == 1:
	
		g_der = classesgamma_der(gam)
		
		probs, probs_der = calcul_probs(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param,der)
		L = np.zeros((1,probs[0][Ancmax].shape[1]))
		L_der = np.zeros((1,probs[0][Ancmax].shape[1]))
		for k in range(4):
			for i,a in ref_root:
				probs[k][Ancmax][:,i] = probs[k][a][:,i]
				probs_der[k][Ancmax][:,i] = probs_der[k][a][:,i]
			Lk = np.sum(probs[k][Ancmax],axis=0)
			L_derk = np.sum(probs_der[k][Ancmax],axis=0)
			L += Lk
			L_der += L_derk
		
		L_der = np.sum(L_der/L)
		L = np.sum(np.log(L/4))
		
		return -L,-L_der
		
	if der == 2:
		probs, probs_der = calcul_probs(probs,tree,BL,roots,coord_root,gam,compos_mean_1,compos_mean_2,eigenvectors_1,eigenvectors_2,code,model_aa,param,der)
		L = np.zeros((1,probs[0][Ancmax].shape[1]))
		L_der = np.zeros((len(param),1,probs[0][Ancmax].shape[1]))
		for k in range(4):
			for i,a in ref_root:
				probs[k][Ancmax][:,i] = probs[k][a][:,i]
				for dd in range(len(param)):
					probs_der[k][Ancmax][dd,:,i] = probs_der[k][a][dd,:,i]
			Lk = np.sum(probs[k][Ancmax],axis=0)
			L_derk = np.sum(probs_der[k][Ancmax],axis=1)
			L += Lk
			for dd in range(len(param)):
				L_der[dd] += L_derk[dd]
			
		L_der = np.sum(L_der/L,axis=0)
		L = np.sum(np.log(L/4))
		
		return -L, -L_der

