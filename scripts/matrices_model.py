import numpy as np
import PCA_module

def mat_gtr(a,b,c,d,e,f,der):
	
	gtr = np.zeros((8,8))
	gtr[:4,:4] = np.array([[0, a, b, c, 0, 0, 0, 0],
						   [a, 0, d, e, 0, 0, 0, 0],
						   [b, d, 0, f, 0, 0, 0, 0],
						   [c, e, f, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, a, b, c],
						   [0, 0, 0, 0, a, 0, d, e],
						   [0, 0, 0, 0, b, d, 0, f],
						   [0, 0, 0, 0, c, e, f, 0]])
	if der == 0:
		return gtr
	elif der == 2:
		gtr_der = [0,0,0,0,0,0]
		gtr_der[0] = np.array([[0, 1, 0, 0, 0, 0, 0, 0],
							   [1, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 1, 0, 0],
							   [0, 0, 0, 0, 1, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0]])
		
		gtr_der[1] = np.array([[0, 0, 1, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [1, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 1, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 1, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0]])

		gtr_der[2] = np.array([[0, 0, 0, 1, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [1, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 1],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 1, 0, 0, 0]])

		gtr_der[3] = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 1, 0, 0, 0, 0, 0],
							   [0, 1, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 1, 0],
							   [0, 0, 0, 0, 0, 1, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0]])

		gtr_der[4] = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 1, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 1, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 1],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 1, 0, 0]])

		gtr_der[5] = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 1, 0, 0, 0, 0],
							   [0, 0, 1, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 0],
							   [0, 0, 0, 0, 0, 0, 0, 1],
							   [0, 0, 0, 0, 0, 0, 1, 0]])
		
		return gtr, gtr_der


def mat_gtr_freq(gtr,freq):
	gtr = freq*gtr
	for i in range(8):
		gtr[i,i] = -sum(gtr[i])
	alp1 = -1/(freq[0]*gtr[0,0]+freq[1]*gtr[1,1]+freq[2]*gtr[2,2]+freq[3]*gtr[3,3])
	gtr[:4,:4] = alp1*gtr[:4,:4]
	alp2 = -1/(freq[4]*gtr[4,4]+freq[5]*gtr[5,5]+freq[6]*gtr[6,6]+freq[7]*gtr[7,7])
	gtr[4:,4:] = alp2*gtr[4:,4:]
	return gtr


def mat_codon(k,w,der):#freq_nt,freq_aa):
	codons = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG',
			  'TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
			  'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG',
			  'CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
			  'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG',
			  'AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
			  'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG',
			  'GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']
	codons_trans = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
					'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X','TGT':'C','TGC':'C','TGA':'X','TGG':'W',
					'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
					'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
					'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

	
	matrix_codons = np.zeros((64,64))
	if der == 2:
		matrix_codons_der1 = np.zeros((64,64))
		matrix_codons_der2 = np.zeros((64,64))
	
	for i,c1 in enumerate(codons):
		for j, c2 in enumerate(codons):
			p = [int(c1[0] != c2[0]), int(c1[1] != c2[1]), int(c1[2] != c2[2])]
			if sum(p) == 1 and codons_trans[c1] != 'X' and codons_trans[c2] != 'X': #potential problem if misalignment with STOP codon
				
				matrix_codons[i,j] = 1
				if codons_trans[c1] != codons_trans[c2]: #non synonymous mutation
					matrix_codons[i,j] = matrix_codons[i,j]*w
					if der == 2:
						matrix_codons_der2[i,j] = 1
					
				if c1[p.index(1)]+c2[p.index(1)] in ['AG','GA','CT','TC']: #transition
					matrix_codons[i,j] = matrix_codons[i,j]*k	
					if der == 2:
						matrix_codons_der1[i,j] = 1
						if codons_trans[c1] != codons_trans[c2]:
							matrix_codons_der1[i,j] = w
							matrix_codons_der2[i,j] = k
	
	if der < 2:
		return matrix_codons
	else:
		return matrix_codons, (matrix_codons_der1, matrix_codons_der2)
	
	
	# 
# 		aa_trans = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
# 			    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'],
# 			    'X': ['TAA', 'TAG', 'TGA'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
# 			    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
# 			    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I': ['ATT', 'ATC', 'ATA'],
# 			    'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'],
# 			    'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
# 			    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
# 			    'G': ['GGT', 'GGC', 'GGA', 'GGG']}
# 	
# 	amino_acid_1 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
# 	
# 	freq_dic = {}
# 	nt = ['A','T','C','G']
# 	for i,c in enumerate(codons):
# 		freq_dic[c] = freq_nt[nt.index(c[0])] * freq_nt[nt.index(c[1])] * freq_nt[4+nt.index(c[2])] #theoretical codon frequency given nucleotide equilibrium at 1st/2nd or 3rd position
# 	freq_dic['TAA'] = 0 #STOP codons frequencies
# 	freq_dic['TAG'] = 0
# 	freq_dic['TGA'] = 0
# 	s = sum(freq_dic.values())
# 	for c in codons:
# 		freq_dic[c] = freq_dic[c]/s
# 	freq_dic2 = {}
# 	for aa in aa_trans:
# 		freq_dic2[aa] = 0
# 		for c in aa_trans[aa]:
# 			freq_dic2[aa] += freq_dic[c]
# 	
# 	for i,c1 in enumerate(codons):
# 		for j, c2 in enumerate(codons):
# 			p = [int(c1[0] != c2[0]), int(c1[1] != c2[1]), int(c1[2] != c2[2])]
# 			if sum(p) == 1 and codons_trans[c1] != 'X' and codons_trans[c2] != 'X': #potential problem if misalignment with STOP codon
# 				
# 				matrix_codons[i,j] = freq_aa[amino_acid_1.index(codons_trans[c2])] * freq_dic[c2] / freq_dic2[codons_trans[c2]] #equilibrium of codon, given expected amino acid frequency and nucleotide frequency (codon bias)
# 				if codons_trans[c1] != codons_trans[c2]: #non synonymous mutation
# 					matrix_codons[i,j] = matrix_codons[i,j]*w
# 				if c1[p.index(1)]+c2[p.index(1)] in ['AG','GA','CT','TC']: #transition
# 					matrix_codons[i,j] = matrix_codons[i,j]*k	
# 				
# 	alp = 0
# 	for i in range(64):
# 		matrix_codons[i,i] = -sum(matrix_codons[i])
# 		alp -= matrix_codons[i,i]
# 	
# 	
# 	return matrix_codons/alp


def mat_codon_freq(mat,freq):
	mat = freq*mat
	for i in range(64):
		mat[i,i] = 0
	
	alp = np.sum(freq*mat.T)
	mat = mat-np.sum(mat,axis=1)*np.identity(64)
	return mat/alp
	


def mat_codon_aa_freq(mat,freq):
	
	mat = freq*mat
	
	codons = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG',
			  'TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
			  'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG',
			  'CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
			  'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG',
			  'AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
			  'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG',
			  'GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']
	codons_trans = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
					'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X','TGT':'C','TGC':'C','TGA':'X','TGG':'W',
					'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
					'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
					'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
	aa_trans = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
			    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'],
			    'X': ['TAA', 'TAG', 'TGA'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
			    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
			    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I': ['ATT', 'ATC', 'ATA'],
			    'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'],
			    'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
			    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
			    'G': ['GGT', 'GGC', 'GGA', 'GGG']}
	
	amino_acid = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
	amino_acid_1 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
	codons_dic = {i:j for i,j in zip(codons,freq)}
	freq_aa = {}
	for aa in amino_acid_1:
		freq_aa[aa] = 0
		for c in aa_trans[aa]:
			freq_aa[aa] += codons_dic[c]
	
	matrix_aa = np.zeros((20,20))
	
	for i,c1 in enumerate(codons):
		for j,c2 in enumerate(codons):
			matrix_aa[amino_acid[codons_trans[c1]],amino_acid[codons_trans[c2]]] += mat[i,j] * freq[j] / freq_aa[codons_trans[c2]]
	
	alp = 0
	for i in range(20):
		matrix_aa[i,i] = 0
	
	alp = np.sum(freq*matrix_aa.T)
	matrix_aa = matrix_aa-np.sum(matrix_aa,axis=1)*np.identity(20)
	
		
	return matrix_aa/alp



def mat_aa(model):
	dic = {}
	
	dicconv = {'model_WAG':'WAG', 'model_JTT':'JTT', 'model_LG':'LG'}
	with open('/Users/pgb/Documents/Data/Roscoff/articles/MDH_thermo_Alvinellidae/analyse_proteome/models.txt','r') as fin:
		lines = fin.readlines()
	for l in lines:
		if 'model' in l:
			m = dicconv[l.strip('\n\r')]
			dic[m] = []
		elif l[0] in ['0','1','2','3','4','5','6','7','8','9'] and l.strip('\n\r')[-1] != ';':
			dic[m].append([float(x) for x in l.strip('\n\r').split(' ') if x != ''])
		
	dicout2 = {}
	for m in dic:
		dicout2[m] = np.zeros((20,20))
		for i,z in enumerate(dic[m]):
			for j,k in enumerate(z):
				dicout2[m][i+1][j] = k
		dicout2[m] = dicout2[m] + dicout2[m].T
	
	return np.array(dicout2[model])


def mat_aa_freq(mat,freq):
	
	mat = freq*mat
	alp = np.sum(freq*mat.T)
	mat = mat-np.sum(mat,axis=1)*np.identity(20)
		
	return mat/alp




def get_freq(compos_1,compos_2,bl,eigenvectors_1,eigenvectors_2,code):
	if code == 'nt' or 'aa':
		compos_1 = PCA_module.get_compositions(compos_1, bl[1:], eigenvectors_1) #equilibrium frequencies on the PCA
	elif code == 'codon' or 'codon_aa':
		if code == 'codon':
			compos_1 = PCA_module.get_compositions(compos_1, bl[1:len(eigenvectors_1)+1], eigenvectors_1) #equilibrium frequencies on the PCA 1 (nt)
			compos_2 = PCA_module.get_compositions(compos_2, bl[1+len(eigenvectors_1):], eigenvectors_2) #equilibrium frequencies on the PCA (aa)
		elif code == 'codon_aa':
			compos_1 = [(1-bl[1])/2, (1-bl[1])/2, bl[1]/2, bl[1]/2, (1-bl[1])/2, (1-bl[1])/2, bl[1]/2, bl[1]/2] #equilibrium frequencies (nt), %GC with %G=%C and %A=%T ATCGATCG
			compos_2 = PCA_module.get_compositions(compos_1, bl[2:], eigenvectors_1) #equilibrium frequencies on the PCA (aa)
		
		#change the 2 nt and aa compositions vector into one describing the expected codons equilibrium frequencies
		
		codons = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG',
			  'TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
			  'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG',
			  'CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
			  'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG',
			  'AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
			  'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG',
			  'GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']
		codons_trans = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
					'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X','TGT':'C','TGC':'C','TGA':'X','TGG':'W',
					'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
					'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
					'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
		aa_trans = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
			    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'],
			    'X': ['TAA', 'TAG', 'TGA'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
			    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
			    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I': ['ATT', 'ATC', 'ATA'],
			    'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'],
			    'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
			    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
			    'G': ['GGT', 'GGC', 'GGA', 'GGG']}
	
		amino_acid_1 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	
		freq_dic = {}
		nt = ['A','T','C','G']
		for i,c in enumerate(codons):
			freq_dic[c] = compos_1[nt.index(c[0])] * freq_nt[nt.index(c[1])] * freq_nt[4+nt.index(c[2])] #theoretical codon frequency given nucleotide equilibrium at 1st/2nd or 3rd position
		freq_dic['TAA'] = 0 #STOP codons frequencies
		freq_dic['TAG'] = 0
		freq_dic['TGA'] = 0
		s = sum(freq_dic.values())
		for c in codons:
			freq_dic[c] = freq_dic[c]/s
		freq_dic2 = {}
		for aa in aa_trans:
			freq_dic2[aa] = 0
			for c in aa_trans[aa]:
				freq_dic2[aa] += freq_dic[c]
	
		compos_3 = []
		for c in codons:
			compos_3.append(compos_2[amino_acid_1.index(codons_trans[c])] * freq_dic[c] / freq_dic2[codons_trans[c]])
	
		compos_1 = [x for x in compos_3]
	
# 	compos_1 = [max(1e-3,x) for x in compos_1]
# 	compos_1 = np.array([x/sum(compos_1) for x in compos_1])
	compos_1 = np.maximum(compos_1,np.ones(len(compos_1))*1e-3)
	compos_1 = compos_1 / np.sum(compos_1)
	return compos_1


def equi_matrix(mat,freq,code):
	if code == 'nt':
		mat = mat_gtr_freq(gtr,freq)
	elif code == 'codon_aa':
		mat = mat_codon_aa_freq(mat,freq)
	elif code == 'codon':
		mat = mat_codon_freq(mat,freq)
	elif code == 'aa':
		mat = mat_aa_freq(mat,freq)
	return mat