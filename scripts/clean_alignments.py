import os

th = 20



dicseq = {}
y = 0
path = './'
for root, dirs, files in os.walk(path, topdown=False):
	for name in files:
		if 'aligned.fa' in name:
			with open(path+name,'r') as fin:
				ll = fin.readlines()
			dicseq[name] = [[],[]]
			for l in ll:
				if l[0] == '>':
					dicseq[name][0].append(l)
					y = 1
				else:
					if y == 0:
						dicseq[name][1][-1]+=l.strip('\n')
					else:
						dicseq[name][1].append(l.strip('\n'))
						y = 0				

for k in dicseq:
	name_file = k.split('_')[:-1]
	name_file = '_'.join(name_file)
	
	for i,l in enumerate(dicseq[k][1]): #enlever portions de moins de th acides aminés entre 2 gaps, 3' et 5' considérés comme gap
		l = list(l+'\n')
		c_begin, c_end = 0,0
		if l[0] != '-':
			c_begin = 0
		for c2,c in enumerate(l[:-1]):
			if l[c2] != '-' and (l[c2+1] == '-' or l[c2+1] == '\n'):
				c_end = c2
			if l[c2] == '-' and l[c2+1] != '-':
				c_begin = c2+1
			
			if c_end-c_begin+1 < th and c_end >= c_begin:
				l[c_begin:c_end+1] = ['-']*(c_end-c_begin+1)
				c_begin = 0
				c_end = 0
		
# 		l = ''.join(l)
		
		dicseq[k][1][i] = l
	
	pos = [] #enlever les positions pour lesquelles on a 1 ou moins acide aminé aligné (non informatif)
	for i,u in enumerate(zip(*dicseq[k][1])):
		uu = ''.join(u)
		if uu.count('-') > len(u)-2:
			pos.append(i)
	
	for i,l in enumerate(dicseq[k][1]):
		l_out = ''.join([l[z] for z in range(len(l)) if z not in pos])
		dicseq[k][1][i] = l_out
	
	
	with open('./'+'_'.join(k.split('_')[:3]),'w') as fout:
		for i,j in zip(dicseq[k][0],dicseq[k][1]):
			fout.write(i)
			fout.write(j)
	

	