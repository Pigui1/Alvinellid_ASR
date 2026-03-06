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
from random import choices


def parse_command_line():
	parser = argparse.ArgumentParser(
		prog="compositional marginal reconstruction",
		description="Compositional phylogeny optimizes phylogenetic trees, taking inton account compositional heterogeneity",
		epilog="--input: MSA, -t/-tree: phylogenetic tree (newick format, rooted), -e/--encode: aa, nt, codon or codon_aa, -a/--axis: number of axis to retain. -c/--constrained: file containing axis and weights; temp for the temperature analysis"
	)

	parser.add_argument(
		'--input', default='', #MSA
	)
	
	parser.add_argument(
		'-e', '--encode',
		default='aa' #aa, nt, codon, codon_aa
	)
	
	parser.add_argument(
		'-m', '--model',
		default='' #optimized parameters
	)
	
	parser.add_argument(
		'-v', '--vectors',
		default='' #eigenvectors
	)
	
	parser.add_argument(
		'-a', '--axis', type=int,
		default=-1 #eigenvectors
	)
	
	parser.add_argument(
		'-o', '--output',
		default='output' #output
	)
	
	parser.add_argument(
		'-asr', '--asr', type=int,
		default=0 #output
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



if __name__ == "__main__":
	args = parse_command_line()
	
	sequences = open_MSA(args['input'])
	
	len_ali = len(list(sequences.values())[0])
	
	#PCA on compositions
	with open(args['vectors'],'r') as fin:
		ll = fin.readlines()
	if args['axis'] == -1:
		args['axis'] = len(ll)-1
	eigenvectors_1 = []
	for l in ll[1:1+args['axis']]:
		eigenvectors_1.append([float(x) for x in l.strip('\n').split('\t')])
	
	
	df_compo = get_composition(sequences,args['encode'])
	compo_means = np.array(pd.DataFrame.mean(df_compo)) #mean sur df pour les compositions de base
	print(df_compo)
	print(len_ali)
	
	with open(args['model'],'r') as fin:
		 ll = fin.readlines()
	gamma = float(ll[1].strip('\n').split('\t')[1])
	coord_root = ll[2].strip('\n').split('\t')[1:]
	coord_root = [float(x) for x in coord_root]
	Ancestors = {}
	BL = {}
	ll = ll[4:]
	for l in ll:
		l = l.strip('\n').split('\t')
		Ancestors[l[0]] = [l[1],l[2],l[3]]
		n = l[0]
		l = [float(x) for x in l[4:]]
		ind = int(len(l)/3)
		BL[n] = [[x for x in l[0:ind]],[x for x in l[ind:2*ind]],[x for x in l[2*ind:3*ind]]]
	
	descendants = []
	for Anc in Ancestors:
		if Anc not in descendants:
			if Ancestors[Anc][0] not in Ancestors:
				descendants.append(Ancestors[Anc][0])
			if Ancestors[Anc][1] not in Ancestors:
				descendants.append(Ancestors[Anc][1])
			if Ancestors[Anc][2] not in Ancestors:
				descendants.append(Ancestors[Anc][2])
			
	
	
# 	seqs_dic = {}
# 	for n in df_compo.index:
# # 		ICI
# # 		weights = np.array(df_compo.loc[n])
# 		seqs_dic[n] = np.array(df_compo.loc[n]).reshape(-1,1)
# 	for l in Ancestors:
# 		seqs_dic[l] = np.zeros((20,len_ali))
	
	seqs_dic = ML_model.get_prob(sequences, Ancestors, args['encode'])
	
	roots = {}
	for n in seqs_dic:
		roots[n] = np.zeros((20,len_ali))
		if n == 'Anc'+str(len(Ancestors)):
			roots[n] = np.ones((20,len_ali))
	ref_root = []
	
	seqs_dic = {0:deepcopy(seqs_dic), 1:deepcopy(seqs_dic), 2:deepcopy(seqs_dic), 3: deepcopy(seqs_dic)} #for each Gamma, 4 categories
	
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
	for an in range(1,len(Ancestors)):
		Anc = 'Anc'+str(an)
		Ancdes2[Anc] = [Anc]
		anc2 = Anc
		while anc2 != 'Anc'+str(len(Ancestors)):
			Ancdes2[Anc].append(Ancestors[anc2][2])
			anc2 = Ancestors[anc2][2]
	
	#Ancestral reconstruction
	
	
	gam = ML_model.classesgamma(gamma)
	seqs_dic = ancestral_reco.calcul_probs(seqs_dic,Ancestors,BL,roots,coord_root,gam,compo_means,None,eigenvectors_1,None,args['encode'],model_aa,None)
	fout = open(args['output']+'_model_'+str(len(coord_root)-1)+'_axes.csv','w')
	residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	fout.write('\t'+'\t'.join(residues)+'\n')
	
	if args['asr'] == 1:
		fout2 = open(args['output']+'_model_asr.csv','w')
	
	
	for Anc in range(len(Ancestors),0,-1):
		Anc_cur = 'Anc'+str(Anc)
		Anc_cut = Ancestors[Anc_cur][0]
		prob_Anc_cut = ancestral_reco.asr(Anc_cut,seqs_dic,Ancestors,BL,descendants,Ancdes,roots,coord_root,gam,compo_means,None,eigenvectors_1,None,args['encode'],model_aa)
		prob_tot = np.zeros(prob_Anc_cut[0].shape)
		for k in range(4):
			prob_Anc_cut[k] = prob_Anc_cut[k]/np.sum(prob_Anc_cut[k],axis=0)
			prob_tot += prob_Anc_cut[k]/4
		
		
		if args['asr'] == 1:
			fout2.write(Anc_cur+'\n')
			for p in prob_tot:
				fout2.write('\t'.join([str(x) for x in p])+'\n')
				
		prob_tot = np.sum(prob_tot, axis = 1)/len_ali
		fout.write(Anc_cur+'\t'+'\t'.join([str(x) for x in prob_tot])+'\n')
	
	fout.close()
	if args['asr'] == 1:
		fout2.close()
	
	