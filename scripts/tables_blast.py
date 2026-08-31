from sys import argv

f1 = argv[1]+'.tabular'
f2 = argv[1] + '_short.csv'

table_short = []
with open(f1,'r') as fin:
	for l in fin.readlines():
		l = l.strip('\n').split('\t')
		table_short.append([l[0],l[1],float(l[3]),float(l[2]),float(l[10])])

table_shorter = []
y = table_short[0][0]
c = [table_short[0]]
for k in table_short:
	if k[0] != y:
		y = k[0]
		toadd = c[0]
		for cc in c:
			if cc[4] < toadd[4]:
				toadd = cc
			elif cc[4] == toadd[4]:
				if cc[2] < toadd[2]:
					toadd = cc
				elif cc[2] == toadd[2]:
					if cc[3] < toadd[3]:
						toadd = cc
		table_shorter.append(toadd)
		c = [k]
	else:
		c.append(k)


toadd = c[0]
for cc in c:
	if cc[4] < toadd[4]:
		toadd = cc
	elif cc[4] == toadd[4]:
		if cc[2] < toadd[2]:
			toadd = cc
		elif cc[2] == toadd[2]:
			if cc[3] < toadd[3]:
				toadd = cc
table_shorter.append(toadd)

diff_seq = []
len_seq = 0
with open(f2,'w') as fout:
	for t in table_shorter:
		fout.write('\t'.join([t[0],t[1],str(int(t[2])),str(t[3]),str(t[4])])+'\n')
		diff_seq.append(t[1])
		len_seq += t[2]

diff_seq = set(diff_seq)
len_seq_avg = len_seq*len(diff_seq)/len(table_shorter)
print('expected match: '+str(len(diff_seq)))
print('expected alignment length: '+str(len_seq_avg))