list_1 = []
with open('sampled_seq_095_1000rep.fa','r') as fin:
	for l in fin:
		if not l.startswith('>'):
			list_1.append(l)


list_2 = []
with open('sampled/sampled_seq_095_1000rep.fa','r') as fin:
	for l in fin:
		if not l.startswith('>'):
			list_2.append(l)


i = 0
for l1 in list_1:
	if l1 not in list_2:
		i+= 1

print(len(list_1),len(list_2))
print(i)