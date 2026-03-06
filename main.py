import argparse
import pandas as pd
import PCA_module
import matrices_model
import ML_model
import ancestral_reco
import tree
import numpy as np
from copy import deepcopy
from scipy.optimize import minimize
import matplotlib.pyplot as plt


def parse_command_line():
	parser = argparse.ArgumentParser(
		prog="compositional phylogeny",
		description="Compositional phylogeny optimizes phylogenetic trees, taking inton account compositional heterogeneity",
		epilog="--input: MSA, -t/-tree: phylogenetic tree (newick format, rooted), -e/--encode: aa, nt, codon or codon_aa, -a/--axis: number of axis to retain. -c/--constrained: file containing axis and weights; temp for the temperature analysis"
	)

	parser.add_argument(
		'--input', default='', #MSA
	)
	
	parser.add_argument(
		'-t', '--tree', #tree
		default=''
	)
	
	parser.add_argument(
		'-e', '--encode',
		default='aa' #aa, nt, codon, codon_aa
	)
	
	parser.add_argument(
		'-a', '--axis', type=int,
		default=-1 #number of axis. 0: BIC evaluation
	)
	
	parser.add_argument(
		'-c', '--constrained',
		default='no' #reads in a file the number/weights of constrained axis. file thermostable
	)
	
	parser.add_argument(
		'-nt', '--threads', type=int,
		default=1 #number of threads, used to parallelized branch length estimation. won't accept more than 4 threads, 2 or 3 advised
	)
	
	parser.add_argument(
		'-o', '--output',
		default='output' #output
	)
	
	parser.add_argument(
		'-asr', '--asr', type=int,
		default=0 #save detailed ancestral reconstruction for each node, each model
	)

	args = parser.parse_args()

	arg_dict = vars(args)
	
	return arg_dict


def open_MSA(msa_file):
	name = ''
	dic_seq = {}
	with open(msa_file, 'r') as f_in:
		for l in f_in:
			if l[0] == '>':
				name = l[1:].strip('\n')
			else:
				dic_seq[name] = l.strip('\n').upper()
	return dic_seq


def get_composition(seqs,res):
	seqs2 = deepcopy(seqs)
	if res == 'aa' or res == 'codon_aa':
		residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	elif res == 'nt' :
		residues = ['A','T','C','G','A3','T3','C3','G3']
	elif res == 'codon':
		residues = ['A','T','C','G','A3','T3','C3','G3']
		residues1 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
	if res == 'codon' or res == 'nt':
		for sp in seqs2:
			i = 2
			seqs2[sp] = list(seqs2[sp])
			while i < len(seqs2[sp]):
				if seqs2[sp][i] != '-':
					seqs2[sp][i] = seqs2[sp][i]+'3'
				i+=3
	dic_compo = {}
	for sp in seqs2:
		dic_compo[sp] = {x:0 for x in residues}
		l_gap = seqs2[sp].count('-')
		len_seq = len(seqs2[sp])
		for r in dic_compo[sp]:
			dic_compo[sp][r] = seqs2[sp].count(r)/(len_seq-l_gap)
	if res == 'nt' or res == 'codon':
		for sp in dic_compo:
			dic_compo[sp]['A'] = dic_compo[sp]['A']*3/2
			dic_compo[sp]['T'] = dic_compo[sp]['T']*3/2
			dic_compo[sp]['C'] = dic_compo[sp]['C']*3/2
			dic_compo[sp]['G'] = dic_compo[sp]['G']*3/2
			dic_compo[sp]['A3'] = dic_compo[sp]['A3']*3
			dic_compo[sp]['T3'] = dic_compo[sp]['T3']*3
			dic_compo[sp]['C3'] = dic_compo[sp]['C3']*3
			dic_compo[sp]['G3'] = dic_compo[sp]['G3']*3
	if res == 'codon':
		codons_trans = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
					    'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X','TGT':'C','TGC':'C','TGA':'X','TGG':'W',
					    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
					    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
					    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
		seqs_prot = {}
		for sp in seqs:
			seqs_prot[sp] = ''
			i = 0
			while i < len(seqs[sp]):
				if '-' in seqs[sp][i:i+3]:
					seqs_prot[sp] += '-'
				else:
					seqs_prot[sp] += codons_trans[seqs[sp][i:i+3]]
				i+=3
		dic_compo_aa = {}
		for sp in seqs2:
			dic_compo_aa[sp] = {x:0 for x in residues1}
			l_gap = seqs_prot[sp].count('-')
			len_seq = len(seqs_prot[sp])
			for res in dic_compo_aa[sp]:
				dic_compo_aa[sp][res] = seqs_prot[sp].count(res)/(len_seq-l_gap)
		return pd.DataFrame.from_dict(dic_compo, orient='index'), pd.DataFrame.from_dict(dic_compo_aa, orient='index')
	else:
		return pd.DataFrame.from_dict(dic_compo, orient='index')


def constrained_PCA(vec_cons,df_compo):

	df_compo_2 = df_compo - pd.DataFrame.mean(df_compo)
	df_compo_3 = df_compo_2 + 0
	
	if vec_cons.ndim == 2: #making orthogonal unit vectors
		for i,_ in enumerate(vec_cons):
			for v0 in vec_cons[:i]:
				coords_proj = vec_cons[i].dot(v0.reshape([-1,1])) #projection of vec_test[i] along v0
				vec_cons[i] = vec_cons[i] - coords_proj*v0 #orthogonal part of vec_test[i]
			vec_cons[i] = vec_cons[i]/np.linalg.norm(vec_cons[i]) #unit vector
	
		for v in vec_cons:
			coords = (df_compo_2.loc[:,:].values).dot(v.reshape([-1,1])) #projection of df_compo along vec_test
			df_compo_3 = df_compo_3 - coords*v #orthogonal residues
			
			cov = (coords.T).dot(coords)/len(df_compo_2.index) #variance (eigenvalue) of component
# 			print(cov)
	
	else: #making unit vector
		vec_cons = vec_cons/np.linalg.norm(vec_cons)
		
		coords = (df_compo_2.loc[:,:].values).dot(vec_cons.reshape([-1,1])) #projection of df_compo along vec_test
		df_compo_3 = df_compo_3 - coords*vec_cons #orthogonal residues
		
		cov = (coords.T).dot(coords)/len(df_compo_2.index) #variance (eigenvalue) of component
# 		print(cov)

	
	return vec_cons, df_compo_3


def writetree(t,BL):
	arbrewrite = {}
	for A in range(1,len(t)+1):
		Anc = 'Anc'+str(A)
		if t[Anc][0] in arbrewrite:
			t1 = arbrewrite[t[Anc][0]]
		else:
			t1 = t[Anc][0]
		if t[Anc][1] in arbrewrite:
			t2 = arbrewrite[t[Anc][1]]
		else:
			t2 = t[Anc][1]
		arbrewrite[Anc] = '('+t1+':'+str(BL[Anc][0][0])+','+t2+':'+str(BL[Anc][1][0])+')'+Anc
		if A == len(t):
			if t[Anc][2] in arbrewrite:
				t3 = arbrewrite[t[Anc][2]]
			else:
				t3 = t[Anc][2]
			arbrewrite[Anc] = '('+arbrewrite[Anc]+':'+str(BL[Anc][2][0]/2)+','+t3+':'+str(BL[Anc][2][0]/2)+');'
	
	return arbrewrite['Anc'+str(len(t))]



if __name__ == "__main__":
	args = parse_command_line()
	
	with open('log_command.txt','w') as f_out:
		f_out.write('log file\n')
	
	sequences = open_MSA(args['input'])
	
	#PCA on compositions
	df_compo = get_composition(sequences,args['encode'])
# 	print(df_compo)
	if args['encode'] == 'codon':
		df_compo_nt = df_compo[0]
		df_compo_aa = df_compo[1]
		eigenvectors_nt = PCA_module.main(df_compo_nt,args['output'])
	# 	print(eigenvectors)
# 		eigenvectors_nt = eigenvectors_nt[:2] #crop sur les eigenvectors retenus
		if args['constrained'] == 'no':
			eigenvectors_aa = PCA_module.main(df_compo_aa,args['output'])
	# 		print(eigenvectors)
# 			eigenvectors_aa = eigenvectors_aa[:2] #crop sur les eigenvectors retenus
		else:
			with open(args['constrained'],'r') as fin:
				vec = []
				for v in fin:
					v = v.strip('\n')
					v = [float(x) for x in v.split('\t')]
					vec.append(v)
			vec = np.array(vec)
			vec_cons, df_compo_aa_2 = constrained_PCA(vec,df_compo_aa) #making orthogonal unit vectors
			eigenvectors_aa = PCA_module.main(df_compo_aa_2,args['output']) #EST-CE QUE MEAN DF_COMPO_2 EST TOUJOURS 0 ?
			eigenvectors_aa = np.vstack((vec_cons,eigenvectors_aa[:-len(vec_cons)]))
	
	
		compo_means_nt = np.array(pd.DataFrame.mean(df_compo_nt)) #mean sur df pour les compositions de base
		comps_nt = PCA_module.get_compositions(compo_means_nt, np.array([[0.4],[0.2]]), eigenvectors_nt)
		compo_means_aa = np.array(pd.DataFrame.mean(df_compo_aa)) #mean sur df pour les compositions de base
		comps_aa = PCA_module.get_compositions(compo_means_aa, np.array([[0.4],[0.2]]), eigenvectors_aa)
	# 	comps = PCA_module.get_compositions(compo_means, np.array(species_pos.loc[1])[:nombre_eigenvectors].reshape([-1,1]), eigenvectors)
	else:
		if args['constrained'] == 'no':
			eigenvectors = PCA_module.main(df_compo,args['output'])
	# 		print(eigenvectors)
# 			eigenvectors = eigenvectors[:2] #crop sur les eigenvectors retenus
		else:
			with open(args['constrained'],'r') as fin:
				ll = fin.readlines()
			vec = []
			for v in ll[1:]:
				v = v.strip('\n')
				v = [float(x) for x in v.split('\t')]
				vec.append(v)
			vec = np.array(vec)
			vec_cons, df_compo_2 = constrained_PCA(vec,df_compo) #making orthogonal unit vectors
			eigenvectors = PCA_module.main(df_compo_2,args['output']) #EST-CE QUE MEAN DF_COMPO_2 EST TOUJOURS 0 ?
			eigenvectors = np.vstack((vec_cons,eigenvectors[:-len(vec_cons)]))
	
	
		compo_means = np.array(pd.DataFrame.mean(df_compo)) #mean sur df pour les compositions de base
# 		comps = PCA_module.get_compositions(compo_means, np.array([[0.4],[0.2]]), eigenvectors)
	# 	comps = PCA_module.get_compositions(compo_means, np.array(species_pos.loc[1])[:nombre_eigenvectors].reshape([-1,1]), eigenvectors)

	
	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	with open(args['output']+'_eigenvectors.txt','w') as fout:
		fout.write('\t'.join(residues)+'\n')
		for e in eigenvectors:
			fout.write('\t'.join([str(x) for x in e])+'\n')
	
	
	
	#OBTENIR LES ARBRES
	
	Ancestors, BL = tree.main(args['tree'])
	descendants = tree.finddesc(tree.gettree(args['tree']))
	
	seqs_dic = ML_model.get_prob(sequences, Ancestors, args['encode'])
	seqs_dic = {0:deepcopy(seqs_dic), 1:deepcopy(seqs_dic), 2:deepcopy(seqs_dic), 3: deepcopy(seqs_dic)} #for each Gamma, 4 categories
	
# 	seqs_dic_der = ML_model.get_prob_der(sequences, Ancestors, args['encode'])
# 	seqs_dic_der = {0:deepcopy(seqs_dic_der), 1:deepcopy(seqs_dic_der), 2:deepcopy(seqs_dic_der), 3: deepcopy(seqs_dic_der)} #for each Gamma, 4 categories
	
# 	seqs_dic_der {0:{x:np.zeros((seqs_dic[0][x].shape)) for x in seqs_dic[0]},1:{x:np.zeros((seqs_dic[0][x].shape)) for x in seqs_dic[0]},2:{x:np.zeros((seqs_dic[0][x].shape)) for x in seqs_dic[0]},3:{x:np.zeros((seqs_dic[0][x].shape)) for x in seqs_dic[0]}}
	
	roots,ref_root = ML_model.get_roots(sequences, Ancestors, args['encode'])
	
	if args['encode'] == 'aa': #special case for the empirical amino acid matrix, which is not optimized
		model_aa = matrices_model.mat_aa('LG')
	
	#incrémentation des eigenvecteurs, minimisation, calcul BIC
	#minimisation paramètres évolution
		#minimisation gamma et BL
	
	
	Ancdes = {}
	for anc in range(1,len(Ancestors)+1):
		Anc = 'Anc'+str(anc)
		Ancdes[Anc] = []
		if Ancestors[Anc][0] not in Ancestors:
			Ancdes[Anc].append(Ancestors[Anc][0])
		else:
			Ancdes[Anc] += [x for x in Ancdes[Ancestors[Anc][0]]]
		if Ancestors[Anc][1] not in Ancestors:
			Ancdes[Anc].append(Ancestors[Anc][1])
		else:
			Ancdes[Anc] += [x for x in Ancdes[Ancestors[Anc][1]]]
		if anc == len(Ancestors):
			if Ancestors[Anc][2] not in Ancestors:
				Ancdes[Anc].append(Ancestors[Anc][2])
			else:
				Ancdes[Anc] += [x for x in Ancdes[Ancestors[Anc][2]]]
	
	Ancdes2 = {}
# 	for anc in range(1,len(Ancestors)+1):
# 		Anc = 'Anc'+str(anc)
# 		Ancdes2[Anc] = []
# 		if Ancestors[Anc][0] not in Ancestors:
# 			Ancdes2[Anc].append(Ancestors[Anc][0])
# 		else:
# 			Ancdes2[Anc] += [Ancestors[Anc][0]]+[x for x in Ancdes2[Ancestors[Anc][0]]]
# 		if Ancestors[Anc][1] not in Ancestors:
# 			Ancdes2[Anc].append(Ancestors[Anc][1])
# 		else:
# 			Ancdes2[Anc] += [Ancestors[Anc][1]]+[x for x in Ancdes2[Ancestors[Anc][1]]]
# 		if anc == len(Ancestors):
# 			if Ancestors[Anc][2] not in Ancestors:
# 				Ancdes2[Anc].append(Ancestors[Anc][2])
# 			else:
# 				Ancdes2[Anc] += [Ancestors[Anc][2]]+[x for x in Ancdes2[Ancestors[Anc][2]]]
	for an in range(1,len(Ancestors)):
		Anc = 'Anc'+str(an)
		Ancdes2[Anc] = [Anc]
		anc2 = Anc
		while anc2 != 'Anc'+str(len(Ancestors)):
			Ancdes2[Anc].append(Ancestors[anc2][2])
			anc2 = Ancestors[anc2][2]

# 	print('OK')
# 	L,L_der = ML_model.calcul_L(seqs_dic,Ancestors,BL,roots,coord_root,gamma,compo_means,None,eigenvectors,None,args['encode'],model_aa,None,ref_root,1)
# 	L = ML_model.calcul_L(seqs_dic,Ancestors,BL,roots,coord_root,gamma+10e-5,compo_means,None,eigenvectors,None,args['encode'],model_aa,None,ref_root,0)
# 	print(L,L_der)
	
# 	ML_model.iter_branch(seqs_dic,Ancestors,BL,descendants,roots,coord_root,gamma,compo_means,None,eigenvectors,None,None,args['encode'],model_aa,ref_root)
	
# 	BL_init = deepcopy(BL)
	
	sample_size = seqs_dic[0]['Anc1'].shape[1]
	
	gamma = 50
	coord_root = [0.0]
	for anc in BL:
		BL[anc][0] = [BL[anc][0]]
		BL[anc][1] = [BL[anc][1]]
		BL[anc][2] = [BL[anc][2]]
	
# 	na = 2
# 	for i in range(1,na):
# 		coord_root.append(0.0)
# 		for anc in BL:
# 			BL[anc][0].append(0.0)
# 			BL[anc][1].append(0.0)
# 			BL[anc][2].append(0.0)
# 	
# 	
# 	for i in range(na,na+1):
	L_vec = []
	AIC_vec = []
	BIC_vec = []
	
	if args['axis'] == -1:
		args['axis'] = 21
	
	for i in range(args['axis']+1):
	
		if args['axis'] == 21 and i >= 3:
			if (AIC_vec[-1] > AIC_vec[-2] and AIC_vec[-2] > AIC_vec[-3]) or i == 20: #we stop the automatic detection when two consecutive AIC are worse
				break
		
		eigenvectors_1 = eigenvectors[:i]
		if i > 0:
			coord_root.append(0.0)
			for anc in BL:
# 				BL[anc][0][0] = BL_init[anc][0]
# 				for j,k in enumerate(BL[anc][0][1:]):
# 					BL[anc][0][j+1] = k
				BL[anc][0].append(0.0)
# 				BL[anc][1][0] = BL_init[anc][1]
# 				for j,k in enumerate(BL[anc][1][1:]):
# 					BL[anc][1][j+1] = k
				BL[anc][1].append(0.0)
# 				BL[anc][2][0] = BL_init[anc][2]
# 				for j,k in enumerate(BL[anc][2][1:]):
# 					BL[anc][2][j+1] = k
				BL[anc][2].append(0.0)
# 		gamma = 50
		print(i)
		
		L = ML_model.calcul_L(gamma,seqs_dic,Ancestors,BL,roots,coord_root,compo_means,None,eigenvectors_1,None,args['encode'],model_aa,None,ref_root,0)
		print(L)
		L_init = 0
		
		while (L_init-L)**2 > 1e-4:
		
			L_init = L
			BL, coord_root = ML_model.iter_branch(L,seqs_dic,Ancestors,BL,descendants,roots,coord_root,gamma,compo_means,None,eigenvectors_1,None,None,args['encode'],model_aa,ref_root,args['threads'])
			
			gamma_res = minimize( ML_model.calcul_L,gamma,(seqs_dic,Ancestors,BL,roots,coord_root,compo_means,None,eigenvectors_1,None,args['encode'],model_aa,None,ref_root,1),method='TNC',jac=True,bounds=[(0.01,50)])
			gamma = gamma_res.x[0]
			
			L = ML_model.calcul_L(gamma,seqs_dic,Ancestors,BL,roots,coord_root,compo_means,None,eigenvectors_1,None,args['encode'],model_aa,None,ref_root,0)
		
		with open(args['output']+'_'+str(i)+'_axis.tre','w') as fout:
			fout.write(writetree(Ancestors,BL))
		
		#Ancestral reconstruction
		
		gam = ML_model.classesgamma(gamma)
		
		seqs_dic = ancestral_reco.calcul_probs(seqs_dic,Ancestors,BL,roots,coord_root,gam,compo_means,None,eigenvectors_1,None,args['encode'],model_aa,None)
		fout = open(args['output']+'_compositions_'+str(i)+'_axes.csv','w')
		residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
		fout.write('\t'+'\t'.join(residues)+'\tnbr_residues\n')
		
		if args['asr'] == 1:
			fout2 = open(args['output']+'_asr_'+str(i)+'_axes.csv','w')
		
		for Anc in range(len(Ancestors),0,-1):
			Anc_cur = 'Anc'+str(Anc)
			Anc_cut = Ancestors[Anc_cur][0]
			prob_Anc_cut = ancestral_reco.asr(Anc_cut,seqs_dic,Ancestors,BL,descendants,Ancdes,roots,coord_root,gam,compo_means,None,eigenvectors_1,None,args['encode'],model_aa)
			prob_tot = np.zeros(prob_Anc_cut[0].shape)
			for k in range(4):
				len_gap = 0
				prob_Anc_cut[k] = prob_Anc_cut[k]/np.sum(prob_Anc_cut[k],axis=0)
				for ind,p in ref_root:
					if p not in Ancdes2[Anc_cut]:
						prob_Anc_cut[k][:,ind] = 0
						len_gap += 1
				prob_tot += prob_Anc_cut[k]/4
			
			if args['asr'] == 1:
				fout2.write(Anc_cur+'\n')
				for p in prob_tot:
					fout2.write('\t'.join([str(x) for x in p])+'\n')
					
			prob_tot = np.sum(prob_tot, axis = 1)
			len_seq_anc = prob_Anc_cut[0].shape[1]-len_gap
			prob_tot = prob_tot/len_seq_anc
			fout.write(Anc_cur+'\t'+'\t'.join([str(x) for x in prob_tot])+'\t'+str(len_seq_anc)+'\n')
		
		fout.close()
		if args['asr'] == 1:
			fout2.close()
			
# 		Anc_cut = 'Anc'+str(len(Ancestors))
# 		prob_tot = np.zeros(seqs_dic[0][Anc_cut].shape)
# 		for k in range(4):
# 			len_gap = 0
# 			seqs_dic[k][Anc_cut] = seqs_dic[k][Anc_cut]/np.sum(seqs_dic[k][Anc_cut],axis=0)
# 			for ind,p in ref_root:
# 				seqs_dic[k][Anc_cut][:,ind] = 0
# 				len_gap += 1
# 			prob_tot += seqs_dic[k][Anc_cut]/4
# 		if args['asr'] == 1:
# 			fout2.write(Anc_cut+'\n')
# 			for p in prob_tot:
# 				fout2.write('\t'.join([str(x) for x in p])+'\n')
# 		prob_tot = np.sum(prob_tot, axis = 1)
# 		len_seq_anc = seqs_dic[0][Anc_cut].shape[1]-len_gap
# 		prob_tot = prob_tot/len_seq_anc
# 		fout.write(Anc_cut+'\t'+'\t'.join([str(x) for x in prob_tot])+'\t'+str(len_seq_anc)+'\n')
# 			#args['asr']
# 		fout.close()
# 		if args['asr'] == 1:
# 			fout2.close()
		
		
		nbr_param = 1 + 2*len(Ancestors)*(1+len(eigenvectors_1)) + (1+len(eigenvectors_1)) + len(eigenvectors_1) #gamma, branches, final branch, root
		BIC_score = nbr_param*np.log(sample_size) + 2*L
		AIC_score = nbr_param*2 + 2*L
# 		print(BL)
# 		print(coord_root)
# 		print(gamma)
		print('BIC: ' + str(BIC_score))
		print('AIC: ' + str(AIC_score))
		L_vec.append(L)
		AIC_vec.append(AIC_score)
		BIC_vec.append(BIC_score)
		
		
		
		with open(args['output']+'_log_command_'+str(len(eigenvectors_1))+'_axis.csv','w') as fout:
			fout.write('Likelihood\t'+str(L)+'\tnbr param\t'+str(nbr_param)+'\tsample size\t'+str(sample_size)+'\n')
			fout.write('gamma\t'+str(gamma)+'\n')
			fout.write('root\t'+'\t'.join([str(x) for x in coord_root])+'\n')
			header = 'branch length\t'+'\t'.join(['coeff. axis '+str(i+1) for i in range(len(eigenvectors_1))])
			fout.write('node\tf1\tf2\tparent\t'+header+'\t'+header+'\t'+header+'\n')
			for Anc in Ancestors:
				fout.write(Anc)
				fout.write('\t')
				fout.write('\t'.join(Ancestors[Anc]))
				fout.write('\t')
				fout.write('\t'.join(['\t'.join([str(k) for k in x]) for x in BL[Anc]]))
				fout.write('\n')
		
		
		if len(eigenvectors_1)>0:
			PCA_module.make_plot(df_compo,eigenvectors_1,None,Ancestors,BL,coord_root,model_aa,gamma,args['encode'],args['output'])

#	import timeit
#	print(timeit.timeit(lambda:m@seqs_dic[nam], number=tt)/tt)
	
	L_vec = -np.array(L_vec) + L_vec[0]
	AIC_vec = -np.array(AIC_vec) + AIC_vec[0]
	BIC_vec = -np.array(BIC_vec) + BIC_vec[0]
	
	plt.clf()
	for x,(i,j,k) in enumerate(zip(L_vec,AIC_vec,BIC_vec)):
		plt.plot([x,x],[0,i],color='k')
		plt.plot([x+0.2,x+0.2],[0,j],color='b')
		plt.plot([x+0.4,x+0.4],[0,k],color='r')
	
	plt.savefig(args['output']+'_L_AIC_BIC.svg')




	
	










