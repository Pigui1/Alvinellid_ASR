import os
from subprocess import run, PIPE
import numpy as np

seq_aligned = {}
seq_not_aligned = {}

for root, dirs, files in os.walk('./', topdown=False):
	for name in files:
		if 'aligned.fa' in name and 'APscaff' in name and 'structure' not in name:
			with open(name,'r') as fin:
				ll = fin.readlines()
			k = 0
			for l in ll:
				if k == 1:
					seq = l
					break
				elif 'Apompejana' in l:
					k = 1
			name2 = name.split('_')[:-1]
			name2 = '_'.join(name2)
			seq_aligned[name2] = seq
		
		if 'aligned.fa' not in name and 'APscaff' in name and 'structure' not in name:
			with open(name,'r') as fin:
				ll = fin.readlines()
			k = 0
			for l in ll:
				if k == 1:
					seq = l
					break
				elif 'Apompejana' in l:
					k = 1
			seq_not_aligned[name] = seq

for k in seq_aligned:
	with open(k+'_structure.fa','w') as fout:
		fout.write('>Apomp_not_aligned\n')
		fout.write(seq_not_aligned[k])
		fout.write('>Apomp_aligned\n')
		fout.write(seq_aligned[k])

	cmd = '/usr/local/bin/mafft --quiet --localpair  --maxiterate 16 --inputorder '+k+'_structure.fa > '+k+'_structure2.fa'
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)

dic_structure = {}
conv_dic = {'E0':'H', 'E1':'S', 'E2':'T', 'B0':'h', 'B1':'s', 'B2':'t'}
with open('689A016D00135532F39D4E72.netsurfp.txt','r') as fin:
	ll = fin.readlines()
for l in ll:
	if l[0] == '#':
		continue
	l = l.split('\t')
	if l[2] not in dic_structure:
		dic_structure[l[2]] = []
	
	dic_structure[l[2]].append(conv_dic[l[0] + str(np.argmax([float(l[-3]),float(l[-2]),float(l[-1])]))])

for k in dic_structure:
	with open(k+'_structure2.fa','r') as fin:
		ll = fin.readlines()
	seq_pomp = ''
	comp = 0
	for l in ll:
		if comp == 2:
			seq_pomp += l.strip()
		if l[0] == '>':
			comp+=1
	dic_structure[k] = [x for i,x in enumerate(dic_structure[k]) if seq_pomp[i]!='-']



for root, dirs, files in os.walk('./', topdown=False):
	for name in files:
		if 'aligned.fa' in name and 'APscaff' in name and 'structure' not in name:
# 			print(name)
			with open(name,'r') as fin:
				ll = fin.readlines()
			file_seq = {}
			for l in ll:
				if l[0] == '>':
					sp = l
				else:
					file_seq[sp] = l.strip()
			seq_ref = file_seq['>Apompejana\n']
			for sp in file_seq:
				file_seq[sp] = ''.join([x for i,x in enumerate(file_seq[sp]) if seq_ref[i] != '-'])
			
			file_seq_H = {}
			file_seq_S = {}
			file_seq_T = {}
			file_seq_h = {}
			file_seq_s = {}
			file_seq_t = {}
			
			name2 = '_'.join(name.split('_')[:-1])
			for sp in file_seq:
				file_seq_H[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='H'])
				file_seq_S[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='S'])
				file_seq_T[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='T'])
				file_seq_h[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='h'])
				file_seq_s[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='s'])
				file_seq_t[sp] = ''.join([file_seq[sp][i] for i,x in enumerate(dic_structure[name2]) if x=='t'])
			
			with open('./HE/'+name2+'_H.fa','w') as fout:
				for sp in file_seq_H:
					fout.write(sp)
					fout.write(file_seq_H[sp]+'\n')
			with open('./SE/'+name2+'_S.fa','w') as fout:
				for sp in file_seq_S:
					fout.write(sp)
					fout.write(file_seq_S[sp]+'\n')
			with open('./TE/'+name2+'_T.fa','w') as fout:
				for sp in file_seq_T:
					fout.write(sp)
					fout.write(file_seq_T[sp]+'\n')
			with open('./hb/'+name2+'_h.fa','w') as fout:
				for sp in file_seq_h:
					fout.write(sp)
					fout.write(file_seq_h[sp]+'\n')
			with open('./sb/'+name2+'_s.fa','w') as fout:
				for sp in file_seq_s:
					fout.write(sp)
					fout.write(file_seq_s[sp]+'\n')
			with open('./tb/'+name2+'_t.fa','w') as fout:
				for sp in file_seq_t:
					fout.write(sp)
					fout.write(file_seq_t[sp]+'\n')

