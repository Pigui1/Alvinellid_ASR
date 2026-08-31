import os

filetow = {'Pectinaria':'','Agunneri':'','Hinvalida':'','Acarldarei':'','Anobothrus':'','Acaudata':'','Apompejana':'','Pfijiensis':'','Pnov':'','Psulfincola':'','Phessleri':'','Pmira':'','Pgrasslei':'','Ppalmiformis':'','Ppandorae':'','Punidentata':'','Mpalmata':'','Skaia':'','Mpalmata':'','Nedwardsi':''}

list_sp = ['Pectinaria','Agunneri','Hinvalida','Acarldarei','Anobothrus','Acaudata','Apompejana','Pfijiensis','Pnov','Psulfincola','Phessleri','Pmira','Pgrasslei','Ppalmiformis','Ppandorae','Punidentata','Mpalmata','Skaia','Nedwardsi']



th = 0

for root, dirs, files in os.walk('./', topdown=False):
	for name in files:
		if 'aligned.fa' in name:
			l_sp = []
			l_seq = []
			l_seq_fij, l_seq_uni, l_seq_hess = '', '', ''
			hess, fij, uni = 0, 0, 0
			with open(name,'r') as fin:
				ll = fin.readlines()
			for l in ll:
				if l[0] == '>':
					if 'hessleri' in l:
						hess = 1
					elif 'fijiensis' in l:
						fij = 1
					elif 'unidentata' in l:
						uni = 1
					else:
						l_sp.append(l.strip('\n')[1:])
				else:
					if hess == 1:
						if len(l.strip('\n')) - l.count('-') > len(l_seq_hess) - l_seq_hess.count('-'):
							l_seq_hess = l.strip('\n')
						hess = 0
					elif fij == 1:
						if len(l.strip('\n')) - l.count('-') > len(l_seq_fij) - l_seq_fij.count('-'):
							l_seq_fij = l.strip('\n')
						fij = 0
					elif uni == 1:
						if len(l.strip('\n')) - l.count('-') > len(l_seq_uni) - l_seq_uni.count('-'):
							l_seq_uni = l.strip('\n')
						uni = 0
					else:
						l_seq.append(l.strip('\n'))
			
			if len(l_seq_fij) > 0:
				l_sp.append('Pfijiensis')
				l_seq.append(l_seq_fij)
			if len(l_seq_uni) > 0:
				l_sp.append('Punidentata')
				l_seq.append(l_seq_uni)
			if len(l_seq_hess) > 0:
				l_sp.append('Phessleri')
				l_seq.append(l_seq_hess)
			
				
			for sp in list_sp:
				if sp not in l_sp:
					l_sp.append(sp)
					l_seq.append('-'*len(l_seq[-1]))
			
			pos = []
			for i,k in enumerate(zip(*l_seq)):
				k2 = ''.join(k)
				if k2.count('-') > th:
					pos.append(i)
			
			
			for i,seq in enumerate(l_seq):
				l_seq[i] = ''.join([x for i,x in enumerate(list(seq)) if i not in pos])
			
			
			
			for sp,seq in zip(l_sp,l_seq):
				filetow[sp] += seq


with open('non_informative_concatenation_0missing.fa','w') as fout:
	for sp in filetow:
		fout.write('>'+sp+'\n')
		fout.write(filetow[sp]+'\n')


