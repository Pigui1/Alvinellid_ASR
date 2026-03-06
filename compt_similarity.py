with open('MDH_expressed.fasta','r') as fin:
	ll = fin.readlines()

seq = []
name = []
for l in ll:
	if l in ll:
		if l[0] == '>':
			name.append(l)
		else:
			seq.append(l)

dicout = {x:{y:0 for y in name} for x in name}
for i,s1 in enumerate(seq):
	for j in range(i+1,len(seq)):
		s2 = seq[j]
		for c1,c2 in zip(s1,s2):
			if c1 == c2:
				dicout[name[i]][name[j]] += 1
				dicout[name[j]][name[i]] += 1

dicout_2 = {}
for s in dicout:
	dicout[s] = [dicout[s][x] for x in dicout[s]]
	dicout_2[s] = min([x for x in dicout[s] if x != 0])
	dicout[s] = sum(dicout[s])

print(dicout)
print(dicout_2)