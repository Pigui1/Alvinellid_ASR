import os

for root, dirs, files in os.walk(".", topdown=False):
	for name in files:
		if 'aligned.fa' in name:
			with open(name,'r') as fin:
				ll = fin.readlines()
			nsp = name.split('_')[0]
			if nsp == 'Nedwardsii':
				nsp = 'Nedwardsi'
			if nsp == 'Terebhyd':
				nsp = 'Skaia'
			nnsp = name.split('_')[1:3]
			nnsp = '_'.join(nsp+nnsp)
			with open(nnsp+'.fa','w') as fout:
				for l in ll:
					if l[0] == '>':
						fout.write(l)
						ide = nnsp.split('_')[0]
						with open('Apomp_'+l[1:].strip('\n')+'_reciprocal_besthits.csv','a') as fout2:
							fout2.write('Apompejana_'+ide+'\t'+l[1:].strip('\n')+'_'+ide+'\t')
					else:
						fout.write(l.replace('-',''))
						with open('Apomp_'+l[1:].strip('\n')+'_reciprocal_besthits.csv','a') as fout2:
							fout2.write(str(len(l.replace('-','').strip('\n')))+'\n')