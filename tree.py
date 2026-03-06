# -- coding: utf-8 --
#!/usr/bin/python


from sys import argv
from copy import deepcopy



def gettree(chemin):
	with open(chemin, 'r') as t:
		tree = t.readline().strip('\n\r')
	return tree
	

def finddepth(tree):
	#checks for the depth of tree nodes, increased if '(' is met, decreased if ',' or ')'
	#returns a list, each character of the tree string is associated to its depth
	depth = 0
	w = 0
	treedepth = []
	for cha in tree:
		if cha == '(':
			depth+=1
		if cha == ',' and w == 1:
			w = 0
			depth-=1
		treedepth.append(depth)
		if cha == ')' and w == 1:
			w = 0
			depth-=1
		if cha == ')':
			if depth == 1:
				treedepth.append(depth)
				break
			else:
				w = 1
	
	maxdepth = max(treedepth)
	return treedepth, maxdepth



def findAnc(tree, treedepth, Ancestors, maxdepth, Anc):
	#creates the value [daughter_node1,daughter_node2,mother_node] in the dictionnary for an entry Anc
	cha = 0
	
	while cha < len(tree):
		Ancname = 'Anc'+str(Anc)
	
		if treedepth[cha] == maxdepth and maxdepth > 1:
			Ancestors[Ancname]=[]
			t = 0
			while treedepth[cha] == maxdepth:
				Ancestors[Ancname].append(tree[cha])
				if tree[cha] == ':':
					pos2d = t #index of the branch length
				treedepth = treedepth[:cha] + treedepth[cha+1:]
				tree = tree[:cha] + tree[cha+1:]
				t+=1
				
				if len(tree) == 0:
					break
			
			bl=''.join(Ancestors[Ancname])[pos2d:]
			tree = tree[:cha]+Ancname+bl+tree[cha:] #the substring (A:,B:)C: is replaced by the new ancestor's name
			treedepth = treedepth[:cha]+[maxdepth-1]*(len(Ancname)+len(bl))+treedepth[cha:]
			
			Ancestors[Ancname] = guessnames(''.join(Ancestors[Ancname])) #transforms the substring into A:xx, B:xx, :xx
			
			for Descendant in Ancestors[Ancname]: #connects the node to its mother node
				D = Descendant.split(':')[0]
				if D in Ancestors:
					Ancestors[D][2] = Ancname+Ancestors[D][2] #goes from ':branch_length' to 'AncX:branch_length'
			
			
			Anc+=1
			
		elif treedepth[cha] == maxdepth and maxdepth == 1: #last node
			Esp = tree.replace('(','').replace(')','').replace(';','').split(',')
			if Ancestors.get(Esp[0].split(':')[0]):
				Ancestors[Esp[0].split(':')[0]][2] = Esp[1].split(':')[0]+':'+str(float(Esp[0].split(':')[1])+float(Esp[1].split(':')[1]))
			if Ancestors.get(Esp[1].split(':')[0]):
				Ancestors[Esp[1].split(':')[0]][2] = Esp[0].split(':')[0]+':'+str(float(Esp[0].split(':')[1])+float(Esp[1].split(':')[1]))
			tree = ''
		
		cha+=1
	
	return tree, Ancestors, treedepth, Anc


def splitAnc(Ancestors):
	#sorts the list A:1, B:2, C:3 in 2 lists A,B,C and 1,2,3 for each Ancestors key
	Anc1 = {}
	AncBL = {}

	for Anc in Ancestors:
		Anc1[Anc] = []
		AncBL[Anc] = []
		for item in Ancestors[Anc]:
			Anc1[Anc].append(item.split(':')[0])
			AncBL[Anc].append(float(item.split(':')[1]))
			
	return Anc1, AncBL


def prodAnc(tree):
	#creates a dictionnary with Ancestors as keys, and [daughter_node1,daughter_node2,mother_node] as value
	Ancestors = {}
	Anc = 1
	treedepth, maxdepth = finddepth(tree)
	while maxdepth >=1:
		tree, Ancestors, treedepth, Anc = findAnc(tree, treedepth, Ancestors, maxdepth, Anc)
		maxdepth-=1
	Ancestors, BL = splitAnc(Ancestors)
	
	return Ancestors, BL


def guessnames(tree):
	#turns the string (A:1,B:2):3 into the list A:1, B:2, :3
	ign = '(),;' #caracteres a ignorer
	
	for i in ign:
		tree = tree.replace(i, '!')
	tree = tree.split('!')
	treeout = []
	
	for g in tree:
		if len(g)>0:
			treeout.append(g)
	
	return treeout


def finddesc(tree):
	#species names in the tree
	desc = guessnames(tree)
	desout = []
	for d in desc:
		if d.split(':')[0] != '':
			desout.append(d.split(':')[0])
	return desout
	

def adapt(tree, Des2, Des):
	#writes a new tree that keeps the initial topology but keeps only the species in the list Des2
	ignored = []
	for D in Des:
		if D not in Des2:
			ignored.append(D)
	
	for i in ignored:
		
		ind = tree.find(i)
		ind2 = ind+len(i)
		while tree[ind2] in [':','0','1','2','3','4','5','6','7','8','9','.','E','-',',']:
			ind2+=1
		if tree[ind-1] == ',':
			ind = ind-1
		tree = tree[:ind]+tree[ind2:]
		
		c = 1
		ind3 = ind
		while c > 0:
			if tree[ind3] == '(':
				c+=1
			if tree[ind3] == ')':
				c-=1
			ind3+=1
		tree = tree[:ind3-1]+tree[ind3:]
		
		ind2 = ind-1
		c = -1
		while c < 0:
			if tree[ind2] == '(':
				c+=1
			if tree[ind2] == ')':
				c-=1
			ind2-=1
		tree = tree[:ind2+1]+tree[ind2+2:]
		
		ind3 -= 2
		
		if tree[ind3] != ';':
			if tree[ind3-1] in ['0','1','2','3','4','5','6','7','8','9'] and tree[ind3+1] in ['0','1','2','3','4','5','6','7','8','9']:
				ind4 = ind3-1
				ind5 = ind3+1
				while tree[ind4] in ['0','1','2','3','4','5','6','7','8','9','.','E','-']:
					ind4-=1
				while tree[ind5] in ['0','1','2','3','4','5','6','7','8','9','.','E','-']:
					ind5+=1
				bl = tree[ind4+1:ind5].split(':')
				bl = str(sum([float(i) for i in bl]))
				tree = tree[:ind4+1]+bl+tree[ind5:]
		else:
			while tree[ind3] != ')':
				ind3-=1
			tree = tree[:ind3+1]+';'

	return tree



def main(chemin):
	tree = gettree(chemin)
	Ancestors, BL = prodAnc(tree)
	return Ancestors, BL

if __name__ == "__main__":
	Ancestors = main(argv[1])

