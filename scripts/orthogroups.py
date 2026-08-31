import os
import numpy as np
from os.path import exists
from os import remove
import matplotlib.pyplot as plt


chromo_dic = {}


#tereb, ampha, apomp, pando, uni, grass, paral

conv1 = {'Skaia':'outgroup1','Nedwardsi':'outgroup1','Pectinaria':'outgroup1','Mpalmata':'outgroup1',
		'Anobothrus':'outgroup3','Acarldarei':'outgroup2','Hinvalida':'outgroup2','Agunneri':'outgroup2',
		'Ppalmiformis':'paralvinella1','Pgrasslei':'paralvinella1',
		'Apompejana':'alvinella',
		'Punidentata1':'unidentata','Punidentata3':'unidentata','Punidentata4':'unidentata',
		'Ppandorae':'pandorae',
		'Pfijiensis1':'paralvinella2','Pfijiensis3':'paralvinella2','Pfijiensis4':'paralvinella2','Pfijiensis5':'paralvinella2',
		'Phessleri1':'paralvinella2','Phessleri3':'paralvinella2','Pmira':'paralvinella2','Psulfincola':'paralvinella2','Pnov':'paralvinella2'}

phylo_conv = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'J':9,'K':10,'L':11}

rmin = 0.8
lmin = 150

# rmin = 0
# lmin = 0

if lmin > 0:
	folder = 'long_or_complete/'
else:
	folder = 'others/'

dicseq = {}
dicl = {}
for root, dirs, files in os.walk(".", topdown=False):
	for name in files:
		if '.csv' in name and 'Acaudata' not in name and 'Pectinaria' not in name:
			with open(name,'r') as fin:
				namesp = name.split('_')[1]
				for l in fin.readlines():
					l = l.strip('\n').split('\t')
					if l[0] not in dicseq:
						dicseq[l[0]] = [l[1]]
						dicl[l[0]] = [int(l[2])]
					else:
						dicseq[l[0]].append(l[1])
						dicl[l[0]].append(int(l[2]))

config_phylos = [0,0,0,0,0,0,0,0,0,0,0,0]
config_lengths = [0,0,0,0,0,0,0,0,0,0,0,0]
missed = 0

with open('../../ref_pompejana/Apomp_vs_REF_short.csv','r') as fin:
	blast_apomp = fin.readlines()
blast_apomp = {x.split('\t')[0]:x.split('\t')[1] for x in blast_apomp}

with open('../../ref_pompejana/gfftot.gff','r') as fin:
	conv_apomp = fin.readlines()
conv_apomp = {x.strip('\n').split('\t')[8]:[x.split('\t')[0],x.split('\t')[3],x.split('\t')[4]] for x in conv_apomp}


pep_seqs = {}
for root, dirs, files in os.walk("../../proteomes/", topdown=False):
	for name in files:
		if 'PEP' in name:
			sp = name.split('_')[0]
			pep_seqs[sp] = {}
			with open('../../proteomes/'+name,'r') as fin:
				for l in fin.readlines():
					if l[0] == '>':
						n = l[1:].strip('\n')
					else:
						pep_seqs[sp][n] = l


for k in dicseq:
	l_part = []
	class_sp = [conv1[x.split('_')[0]] for x in dicseq[k]]
	cons = 0
	if 'outgroup1' not in class_sp and 'outgroup3' in class_sp:
		l_anob = 0
		l_amphar = []
		for u,v in zip(class_sp,dicl[k]):
			if u == 'outgroup3':
				l_anob = v
			if u == 'outgroup2':
				l_amphar.append(v)
		cons = 1
		if len(l_amphar) >= 2:
			if l_anob <= sorted(l_amphar)[-2]: #Anobothrus ne sera pas conservé
				cons = 0
	if cons == 0:
		for u,v in zip(dicseq[k],dicl[k]):
			if 'Anobothrus' not in u:
				l_part.append(v)
	else:	
		l_part = [x for x in dicl[k]]
	l_max = max(l_part) #longueur maximum alignée, Anobothrus non compté car souvent retiré
	
	sp = []
	ll = []
	for u,v in zip(dicseq[k],dicl[k]):
		if v/l_max > rmin or v > lmin: #longueur similaires des séquences ou séquence suffisamment longue
			sp.append(u)
			ll.append(v)
	dicseq[k] = [x for x in sp]
	dicl[k] = [x for x in ll]
	
	outgroup1, outgroup1_l, outgroup1_n = [0,0], [0,0], [0,0] #terebbelidae
	outgroup2, outgroup2_l, outgroup2_n = [0,0], [0,0], [0,0] #ampharetidae
	outgroup3, outgroup3_l, outgroup3_n = 0, 0, 0
	paralvinella1, paralvinella1_l, paralvinella1_n = 0, 0, 0 #palmiformis ou grasslei
	paralvinella2, paralvinella2_l, paralvinella2_n = 0, 0, 0 #autre paralvinella
	pandorae, pandorae_l, pandorae_n = 0, 0, 0
	unidentata, unidentata_l, unidentata_n = 0, 0, 0
	
	for u, v in zip(dicseq[k],dicl[k]):
		n = u.split('_')[0]
		if conv1[n] == 'outgroup1':
			if v > min(outgroup1_l):
				i = outgroup1_l.index(min(outgroup1_l))
				outgroup1[i] = u
				outgroup1_l[i] = v
				outgroup1_n[i] = 1
		elif conv1[n] == 'outgroup2':
			if v > min(outgroup2_l):
				i = outgroup2_l.index(min(outgroup2_l))
				outgroup2[i] = u
				outgroup2_l[i] = v
				outgroup2_n[i] = 1
		elif conv1[n] == 'outgroup3':
			if v >outgroup3_l:
				outgroup3 = u
				outgroup3_l = v
				outgroup3_n = 1
		elif conv1[n] == 'paralvinella1':
			if v > paralvinella1_l:
				paralvinella1 = u
				paralvinella1_l = v
				paralvinella1_n = 1
		elif conv1[n] == 'paralvinella2':
			if v > paralvinella2_l:
				paralvinella2 = u
				paralvinella2_l = v
				paralvinella2_n = 1
		elif conv1[n] == 'pandorae':
			if v > pandorae_l:
				pandorae = u
				pandorae_l = v
				pandorae_n = 1
		elif conv1[n] == 'unidentata':
			if v > unidentata_l:
				unidentata = u
				unidentata_l = v
				unidentata_n = 1
	
	l_tot = outgroup1_l+outgroup2_l+[paralvinella1_l,paralvinella2_l,pandorae_l,unidentata_l]
	l_tot = [x for x in l_tot if x != 0]
	
	
	tereb = sum(outgroup1_n)
	amphar = sum(outgroup2_n)
	anob = outgroup3_n
	para1 = paralvinella1_n
	para2 = paralvinella2_n
	pando = pandorae_n
	uni = unidentata_n
	
	phylo = 0
	
	if pando == 1 and uni == 1 and para1+para2 >=1: #condition minimale
		if tereb == 0: #si il n'y a pas de Terebbelid, alors Anobothrus peut être considéré comme une espèce outgroup
			if outgroup3_l > min(outgroup2_l):
				i = outgroup2_l.index(min(outgroup2_l))
				outgroup2[i] = outgroup3
				outgroup2_l[i] = outgroup3_l
				l_tot = outgroup1_l+outgroup2_l+[paralvinella1_l,paralvinella2_l,pandorae_l,unidentata_l]
				l_tot = [x for x in l_tot if x != 0]
				outgroup2_n[i] = 1
				amphar = sum(outgroup2_n)
		if tereb == 2 and amphar == 2 and para1+para2 == 1:
			phylo = 'F'
		if tereb == 2 and amphar == 2 and para1+para2 == 2:
			phylo = 'A'
		if tereb == 1 and amphar == 2 and para1+para2 == 1:
			phylo = 'G'
		if tereb == 1 and amphar == 2 and para1+para2 == 2:
			phylo = 'B'
		if tereb == 2 and amphar == 1 and para1+para2 == 1:
			phylo = 'H'
		if tereb == 2 and amphar == 1 and para1+para2 == 2:
			phylo = 'C'
		if tereb+amphar == 2 and para1+para2 == 1:
			phylo = 'I'
		if tereb+amphar == 2 and para1+para2 == 2:
			phylo = 'D'
		if tereb+amphar == 1 and para1+para2 == 1:
			phylo = 'J'
		if tereb+amphar == 1 and para1+para2 == 2:
			phylo = 'E'
		if tereb+amphar == 0 and para1+para2 == 1:
			phylo = 'L'
		if tereb+amphar == 0 and para1+para2 == 2:
			phylo = 'K'
		
		
		if k in blast_apomp:
			
			config_phylos[phylo_conv[phylo]] +=1
			config_lengths[phylo_conv[phylo]] += min(l_tot)
			
			print(outgroup1, outgroup2,paralvinella1,paralvinella2,pandorae,unidentata,k,'-'.join(conv_apomp[blast_apomp[k]]))
			
			with open('../../gene_alignments/'+folder+'-'.join(conv_apomp[blast_apomp[k]])+'_'+phylo+'.fa','w') as fout: #nom pour l'alignement
				fout.write('>Apompejana\n')
				fout.write(pep_seqs[k.split('_')[0]][k])
				
				fout.write('>Ppandorae\n')
				fout.write(pep_seqs[pandorae.split('_')[0]][pandorae])
				
				fout.write('>Punidentata\n')
				fout.write(pep_seqs[unidentata.split('_')[0]][unidentata])
				
				if para1+para2 == 2:
					fout.write('>Pparalvinella1\n')
					fout.write(pep_seqs[paralvinella1.split('_')[0]][paralvinella1])
					fout.write('>Pparalvinella2\n')
					fout.write(pep_seqs[paralvinella2.split('_')[0]][paralvinella2])
				else:
					fout.write('>Pparalvinella\n')
					if para1 == 1:
						fout.write(pep_seqs[paralvinella1.split('_')[0]][paralvinella1])
					elif para2 == 1:
						fout.write(pep_seqs[paralvinella2.split('_')[0]][paralvinella2])
				if tereb+amphar >= 3:
					if tereb == 2:
						fout.write('>Outgroup1\n')
						fout.write(pep_seqs[outgroup1[0].split('_')[0]][outgroup1[0]])
						fout.write('>Outgroup2\n')
						fout.write(pep_seqs[outgroup1[1].split('_')[0]][outgroup1[1]])
						
						fout.write('>Outgroup3\n')
						fout.write(pep_seqs[outgroup2[0].split('_')[0]][outgroup2[0]])
						if amphar == 2:
							fout.write('>Outgroup4\n')
							fout.write(pep_seqs[outgroup2[1].split('_')[0]][outgroup2[1]])
					else:
						fout.write('>Outgroup1\n')
						fout.write(pep_seqs[outgroup2[0].split('_')[0]][outgroup2[0]])
						fout.write('>Outgroup2\n')
						fout.write(pep_seqs[outgroup2[1].split('_')[0]][outgroup2[1]])
						
						fout.write('>Outgroup3\n')
						fout.write(pep_seqs[outgroup1[0].split('_')[0]][outgroup1[0]])
				else:
					nbr_outgroup = 0
					if tereb > 0:
						nbr_outgroup += 1
						fout.write('>Outgroup'+str(nbr_outgroup)+'\n')
						fout.write(pep_seqs[outgroup1[0].split('_')[0]][outgroup1[0]])
						if tereb > 1:
							nbr_outgroup += 1
							fout.write('>Outgroup'+str(nbr_outgroup)+'\n')
							fout.write(pep_seqs[outgroup1[1].split('_')[0]][outgroup1[1]])
					if amphar > 0:
						nbr_outgroup += 1
						fout.write('>Outgroup'+str(nbr_outgroup)+'\n')
						fout.write(pep_seqs[outgroup2[0].split('_')[0]][outgroup2[0]])
						if amphar > 1:
							nbr_outgroup += 1
							fout.write('>Outgroup'+str(nbr_outgroup)+'\n')
							fout.write(pep_seqs[outgroup2[1].split('_')[0]][outgroup2[1]])
				
			
			if folder == 'long_or_complete/':
				for conf in ['A','B','C','D','E','F','G','H','I','J','K','L']:
					if exists('../../gene_alignments/others/'+'-'.join(conv_apomp[blast_apomp[k]])+'_'+conf+'.fa'):
						os.remove('../../gene_alignments/others/'+'-'.join(conv_apomp[blast_apomp[k]])+'_'+conf+'.fa')
			
			if conv_apomp[blast_apomp[k]][0] not in chromo_dic:
				chromo_dic[conv_apomp[blast_apomp[k]][0]] = [conv_apomp[blast_apomp[k]][1]]
			else:
				chromo_dic[conv_apomp[blast_apomp[k]][0]].append(conv_apomp[blast_apomp[k]][1])
		else:
			missed += 1


print(config_phylos, sum(config_phylos))
print(config_lengths, sum(config_lengths))
# print(missed)


i = 0
for k in chromo_dic:
	i+=1
	y = [float(x) for x in chromo_dic[k]]
	plt.scatter([i]*len(chromo_dic[k]),y)
plt.show()

