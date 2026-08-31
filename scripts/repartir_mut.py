with open('individual_list.txt','r') as fin:
	ll = fin.readlines()

i = len(ll)/16
p = 0

for j,l in enumerate(ll):
	if j+1 > i*p:
		p+=1
	with open('seq'+str(p)+'/individual_list_095_1000rep.txt','a') as fout:
		fout.write(l)