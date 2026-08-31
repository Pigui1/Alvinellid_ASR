from sys import argv

f1 = argv[1]+'_vs_'+argv[2]+'_short.csv'
f2 = argv[2]+'_vs_'+argv[1]+'_short.csv'


with open(f1,'r') as fin:
	ll1 = fin.readlines()

with open(f2,'r') as fin:
	ll2 = fin.readlines()


dic1 = {}
for l in ll1:
	l = l.split('\t')
	dic1[l[0]] = [l[1],int(l[2])]

dic2 = {}
for l in ll2:
	l = l.split('\t')
	dic2[l[0]] = [l[1],int(l[2])]

towrite = []
len_tot = 0
len_short = 0
nbr_genes_short = 0
nbr_genes = 0
for n in dic1:
	if dic1[n][0] in dic2:
		if n == dic2[dic1[n][0]][0]:
			towrite.append(n+'\t'+dic1[n][0]+'\t'+str(dic1[n][1])+'\n')
			nbr_genes += 1
			len_tot += dic1[n][1]
			if dic1[n][1] < 150:
				len_short += dic1[n][1]
				nbr_genes_short += 1



with open(argv[1]+'_'+argv[2]+'_reciprocal_besthits.csv','w') as fout:
	for n in towrite:
		fout.write(n)

print('number of reciprocal hits and length: '+str(nbr_genes) + ', '+str(len_tot))
print('number of reciprocal hits and length under 150 amino acids: '+str(nbr_genes_short) + ', '+str(len_short))