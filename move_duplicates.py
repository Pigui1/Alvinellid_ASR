import os

def starting(v):
	return int(v.split('_')[1].split('-')[1])
def stopping(v):
	return int(v.split('_')[1].split('-')[2])
def chroming(v):
	return v.split('_')[0]+'_'+n.split('_')[1].split('-')[0]
def phyling(v):
	return v.split('_')[-1].split('.')[0]
	

folder1 = 'long_or_complete/'
folder2 = 'others/'
folder3 = 'oldgenes/'

convphy = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8,'I':9,'J':10,'K':11,'L':12}

names1 = []
path = './'+folder1
for root, dirs, files in os.walk(path, topdown=False):
	for name in files:
		if '.fa' in name and '_aligned' not in name and 'iqtree' not in name:
			names1.append(name)

dic_ref1 = {}
for n in names1:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	phylo = phyling(n)
	if chromo not in dic_ref1:
		dic_ref1[chromo] = []
	dic_ref1[chromo].append([start,stop,phylo])

names2 = []
path = './'+folder2
for root, dirs, files in os.walk(path, topdown=False):
	for name in files:
		if '.fa' in name and '_aligned' not in name and 'iqtree' not in name:
			names2.append(name)

dic_ref2 = {}
for n in names2:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	phylo = phyling(n)
	if chromo not in dic_ref2:
		dic_ref2[chromo] = []
	dic_ref2[chromo].append([start,stop,phylo])

names3 = []
path = './'+folder3
for root, dirs, files in os.walk(path, topdown=False):
	for name in files:
		if '.fa' in name and '_aligned' not in name and 'iqtree' not in name:
			names3.append(name)

dic_ref3 = {}
for n in names3:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	phylo = phyling(n)
	if chromo not in dic_ref3:
		dic_ref3[chromo] = []
	dic_ref3[chromo].append([start,stop,phylo])


# comparaison long vs long
# doublons1 = 0
# dic_compt = {}
# for n in names1:
# 	chromo = chroming(n)
# 	start = starting(n)
# 	stop = stopping(n)
# 	phylo = phyling(n)
# 	if chromo+'_start'+str(start) not in dic_compt:
# 		dic_compt[chromo+'_start'+str(start)] = 0
# 	dic_compt[chromo+'_start'+str(start)] += 1
# 	if chromo+'_stop'+str(stop) not in dic_compt:
# 		dic_compt[chromo+'_stop'+str(stop)] = 0
# 	dic_compt[chromo+'_stop'+str(stop)] += 1
# 	for k in dic_ref1[chromo]:
# 		if ((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1) or (stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1)) and not (start == k[0] and stop == k[1] and phylo == k[2]):
# 			os.rename('./'+folder1+n, './'+folder1+'doublons/'+n)
# 			doublons1 += 1
# 			if convphy[phylo] > convphy[k[2]]:
# 				print(n,chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa')
# 				os.rename('./'+folder1+n, './doublonslong/'+n)
# 			else:
# 				print(chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa',n)
# 				os.rename('./'+folder1+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa', './doublonslong/'+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa')
# 				print('ok')
# 			print((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1),(stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1))
# print(doublons1)
# 
# exit()
# 
# comparaison short vs short
# doublons1 = 0
# dic_compt = {}
# for n in names2:
# 	chromo = chroming(n)
# 	start = starting(n)
# 	stop = stopping(n)
# 	phylo = phyling(n)
# 	if chromo+'_start'+str(start) not in dic_compt:
# 		dic_compt[chromo+'_start'+str(start)] = 0
# 	dic_compt[chromo+'_start'+str(start)] += 1
# 	if chromo+'_stop'+str(stop) not in dic_compt:
# 		dic_compt[chromo+'_stop'+str(stop)] = 0
# 	dic_compt[chromo+'_stop'+str(stop)] += 1
# 	for k in dic_ref2[chromo]:
# 		if ((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1) or (stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1)) and not (start == k[0] and stop == k[1] and phylo == k[2]):
# 			os.rename('./'+folder2+n, './'+folder2+'doublons/'+n)
# 			doublons1 += 1
# 			if convphy[phylo] > convphy[k[2]]:
# 				print(n,chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa')
# 				os.rename('./'+folder2+n, './doublonsothers/'+n)
# 			else:
# 				print(chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa',n)
# 				os.rename('./'+folder2+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa', './doublonsothers/'+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa')
# 			print((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1),(stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1))
# print(doublons1)
# 
# comparaison old vs old
# doublons1 = 0
# dic_compt = {}
# for n in names3:
# 	chromo = chroming(n)
# 	start = starting(n)
# 	stop = stopping(n)
# 	phylo = phyling(n)
# 	if chromo+'_start'+str(start) not in dic_compt:
# 		dic_compt[chromo+'_start'+str(start)] = 0
# 	dic_compt[chromo+'_start'+str(start)] += 1
# 	if chromo+'_stop'+str(stop) not in dic_compt:
# 		dic_compt[chromo+'_stop'+str(stop)] = 0
# 	dic_compt[chromo+'_stop'+str(stop)] += 1
# 	for k in dic_ref3[chromo]:
# 		if ((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1) or (stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1)) and not (start == k[0] and stop == k[1] and phylo == k[2]):
# 			os.rename('./'+folder3+n, './'+folder3+'doublons/'+n)
# 			doublons1 += 1
# 			if convphy[phylo] > convphy[k[2]]:
# 				print(n,chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa')
# 				os.rename('./'+folder3+n, './doublonsold/'+n)
# 			else:
# 				print(chromo+'-'+str(start)+'-'+str(stop)+'_'+k[2]+'.fa',n)
# 				os.rename('./'+folder3+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa', './doublonsold/'+chromo+'-'+str(k[0])+'-'+str(k[1])+'_'+k[2]+'.fa')
# 			print((start == k[0] and dic_compt[chromo+'_start'+str(start)]>1),(stop == k[1] and dic_compt[chromo+'_stop'+str(stop)]>1))
# print(doublons1)
# 
# 
# exit()


#comparaison long vs short
doublons2 = 0
for n in names2:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	for k in dic_ref1[chromo]:
		if (start >= k[0] and start <= k[1]) or (stop >= k[0] and stop <= k[1]):
			os.rename('./'+folder2+n, './doublonsothers/'+n)
			print(n)
			doublons2 += 1
print(doublons2)
	
	
#comparaison long vs old
doublons3 = 0
for n in names3:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	for k in dic_ref1[chromo]:
		if (start >= k[0] and start <= k[1]) or (stop >= k[0] and stop <= k[1]):
			os.rename('./'+folder3+n, './doublonsold/'+n)
			print(n)
			doublons3 += 1
print(doublons3)


#comparaison short vs old
doublons4 = 0
for n in names3:
	chromo = chroming(n)
	start = starting(n)
	stop = stopping(n)
	if chromo in dic_ref2:
		for k in dic_ref2[chromo]:
			if (start >= k[0] and start <= k[1]) or (stop >= k[0] and stop <= k[1]):
				os.rename('./'+folder3+n, './doublonsold/'+n)
				doublons4 += 1
				print(n)
print(doublons4)










