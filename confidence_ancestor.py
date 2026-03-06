from subprocess import run, PIPE
from sys import argv
import numpy as np
from scipy.optimize import curve_fit,minimize
import matplotlib.pyplot as plt
from numpy.random import choice

#python confidence_ancestor modèle_swiss protéine_référence séquences_fasta fichier_thermodynamique

#python file_marginal_probabilities file_alignement_used_for_ASR name_ancestor_in_maginal_probabilities name_ancestor_in_alignement clean_alignment_file_with_ancestral_sequences modèle_swiss

#modéliser swiss modeller Anc1

#foldx repair jusqu'à convergence
def Repair_model(prot,th):
	print('Optimize model')
	cmd = './foldxMacC11/foldx_20241231 --command=RepairPDB  --pdb='+ prot
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
	std = std.stdout.split('\n')
	for i,l in enumerate(std):
		if 'End of Repair' in l:
			score = float(std[i-2].split('\t')[-1])
			print(score)

	score2 = score-10
	while (score2-score)**2 > th:
		score2 = score
		cmd = './foldxMacC11/foldx_20241231 --command=RepairPDB  --pdb='+prot.split('.')[0]+'_Repair.pdb'
		std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
		std = std.stdout.split('\n')
		for i,l in enumerate(std):
			if 'End of Repair' in l:
				score = float(std[i-2].split('\t')[-1])
				print(score)
		std = run('mv '+prot.split('.')[0]+'_Repair_Repair.pdb '+prot.split('.')[0]+'_Repair.pdb', stdout=PIPE, stderr=PIPE, shell=True, text=True)
	
	cmd = './foldxMacC11/foldx_20241231 --command=Stability --pdb='+prot.split('.')[0]+'_Repair.pdb'
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
	std = std.stdout.split('\n')
	for i,l in enumerate(std):
		if 'FINISHING STABILITY ANALYSIS OPTION' in l:
			score = float(std[i-2].split('\t')[-1])
	
	print('OK')
	return score

#mutations pour les autres séquences

def list_mut(prot, seqin):
	
	with open(seqin,'r') as fin:
		ll = fin.readlines()

	sequences = []
	name_seq = []
	w = 1
	for l in ll:
		if l[0] != '>' and w == 1:
			sequences.append(l.strip('\n'))
		elif l[0] != '>' and w == 0:
			wt = l.strip('\n')
		else:
			if '>'+prot+'\n' == l:
				w = 0
			else:
				w = 1
				name_seq.append(l[1:].strip('\n'))
	
	fout = open('individual_list.txt','a')
	for s in sequences:
		i = 1
		top = ''
		while i <= len(wt):
			if wt[i-1] != s[i-1] and wt[i-1] != '-' and s[i-1] != '-' and i <= 137: #le modèle de swissmodel n'a pas intégré les 4 derniers acides aminés dans la conformation 3D
				if len(top)>0:
					top += ','
# 				top += wt[i-1]+'A'+str(i)+s[i-1]+','+wt[i-1]+'B'+str(i)+s[i-1]
				top += wt[i-1]+'A'+str(i)+s[i-1]
			i+=1
		if len(top) == 0:
			top = wt[0]+'A'+str(1)+wt[0]
		top += ';\n'
		fout.write(top)
	fout.close()
	
	return name_seq



#foldx buildmodel

def build_model(modele_init,prots,rename):
	print('Building models')
	cmd = './foldxMacC11/foldx_20241231 --command=BuildModel --pdb='+modele_init.split('.')[0]+'_Repair.pdb --mutant-file=individual_list.txt'
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
# 	print(std.stdout)
	print('OK\nComputing stabilities')
	k = 0
	scores = []
	for p in prots:
		k+=1
		if rename == 1:
			std = run('mv '+modele_init.split('.')[0]+'_Repair_'+str(k)+'.pdb '+modele_init.split('.')[0]+'_Repair_'+p+'.pdb', stdout=PIPE, stderr=PIPE, shell=True, text=True)
			cmd = './foldxMacC11/foldx_20241231 --command=Stability --pdb='+modele_init.split('.')[0]+'_Repair_'+p+'.pdb'
		else:
			cmd = './foldxMacC11/foldx_20241231 --command=Stability --pdb='+modele_init.split('.')[0]+'_Repair_'+str(k)+'.pdb'
		std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
		std = std.stdout.split('\n')
		for i,l in enumerate(std):
			if 'FINISHING STABILITY ANALYSIS OPTION' in l:
				score = float(std[i-2].split('\t')[-1])
		scores.append(score)
	print('OK')
	return scores
	

#corrélation buildmodel - mesures, optimiser température


def cor_r2(temp,xdata,dH,dCp,Tm,prot,filtre):
	xdata = np.array(xdata)
	ydata = []
	temp = temp[0]
	for k in prot:
		ydata.append(gibbs_helmotz(temp,dH[k],dCp[k],Tm[k]))
	ydata = np.array(ydata) - ydata[-1]
	xd = []
	yd = []
	for i,k in enumerate(prot):
		if k not in filtre:
			xd.append(xdata[i])
			yd.append(ydata[i])
	xd = np.array(xd)
	yd = np.array(yd)
	popt, pcov = curve_fit(lambda x,a:a*x, xd, yd)
	residuals = yd - popt[0]*xd
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((yd-np.mean(yd))**2)
	r_squared = 1 - (ss_res / ss_tot)
	return 1-r_squared


def gibbs_helmotz(T,dH,dCp,Tm):
	Tm = Tm+273.15
	T = T+273.15
	return (dH*(1-T/Tm) - dCp*(Tm-T+T*np.log(T/Tm)))/4186
	

def cor_r2_2(temp,xdata,dH,dCp,Tm,prot,filtre):
	xdata = np.array(xdata)
	ydata = []
	for k in prot:
		ydata.append(gibbs_helmotz(temp,dH[k],dCp[k],Tm[k]))
	ydata = np.array(ydata) - ydata[-1]
	
	xd = []
	yd = []
	pd = []
	for i,k in enumerate(prot):
		if k not in filtre: #filtrer outliers
			xd.append(xdata[i])
			yd.append(ydata[i])
			pd.append(k)
	xd = np.array(xd)
	yd = np.array(yd)
	popt, pcov = curve_fit(lambda x,a:a*x, xd, yd) #forcée ici par 0,0
	residuals = yd - popt[0]*xd
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((yd-np.mean(yd))**2)
	r_squared = 1 - (ss_res / ss_tot)
	return r_squared, popt[0], xdata, ydata


def optim_corr(scores,file_thermo,prot_names):
	dicnames = {'Pu':'Punidentata','Pk':'Ppandorae','Pg':'Pgrasslei','Pp':'Ppalmiformis','Ap':'Apompejana','Pf':'Pfijiensis','Ps':'Psulfincola','Pm':'Pmira','Anc1':'Anc1','Anc2':'Anc2','Anc3':'Anc3','Anc4':'Anc4','Anc5':'Anc5','Anc6':'Anc6','Anc1-T9':'Anc1-T9','Anc2-T9':'Anc2-T9','Anc3-T9':'Anc3-T9','Anc4-T9':'Anc4-T9','Anc5-T9':'Anc5-T9','Anc6-T9':'Anc6-T9'}
	dicnames = {'Pu1':'Punidentata1','Pk1':'Ppandorae1','Pg1':'Pgrasslei1','Pg2':'Pgrasslei2','Pp1':'Ppalmiformis1','Pp2':'Ppalmiformis2','Ap1':'Apompejana1','Ap2':'Apompejana2','Pf':'Pfijiensis','Ps2':'Psulfincola2','Pm1':'Pmira1','Anc1':'Anc1','Anc2':'Anc2','Anc3':'Anc3','Anc4':'Anc4','Anc5':'Anc5','Anc6':'Anc6','Anc1-T9':'Anc1-T9','Anc2-T9':'Anc2-T9','Anc3-T9':'Anc3-T9','Anc4-T9':'Anc4-T9','Anc5-T9':'Anc5-T9','Anc6-T9':'Anc6-T9'}
	filtre = ['Psulfincola','Pgrasslei1','Anc4','Ppalmiformis1','Anc2']
	filtre = ['Anc4']
	
	with open(file_thermo,'r') as fin:
		ll = fin.readlines()
	dH = {}
	dCp = {}
	Tm = {}
	for l in ll:
		l = l.strip('\n').split('\t')
# 		if l[0] != 'Pm' and l[0] != 'Anc2' and l[0] != 'Anc4' and l[0] != 'Anc2-T9' and l[0] != 'Anc3-T9' and l[0] != 'Anc4-T9' and l[0] != 'Anc5-T9':
		dH[dicnames[l[0]]] = float(l[1])*1000
		dCp[dicnames[l[0]]] = float(l[2])*1000
		Tm[dicnames[l[0]]] = float(l[3])
	
	tempopt = minimize(cor_r2, 70, args=(scores,dH,dS,prot_names,filtre), method='TNC', bounds=[(10,100)])
	print('Optimal temperature: '+str(tempopt.x[0]))
	
	tempopt.x[0] = 60
	r_square, a, xd, yd = cor_r2_2(tempopt.x[0],scores,dH,dCp,Tm,prot_names,filtre)
	plt.scatter(xd,yd)
	for p,x,y in zip(prot_names,xd,yd):
		plt.annotate(p,(x,y))
	plt.plot([np.min(xd),np.max(xd)],[a*np.min(xd),a*np.max(xd)])
	print(r_square,a)
	plt.show()
	
	return a,yd
	
	
	
	
#sampling 100 séquences Anc1 (prob > 0.9 considérées certaines)

def sampling_Anc(file_prob,file_alignement,ancestor_file_prob,alignement_ancestor,sampling_number,th):
	
	with open(file_alignement,'r') as fin:
		ll = fin.readlines()
	a = 0
	seq = ''
	for l in ll:
		if a == 1:
			if l[0] == '>':
				break
			seq += l.strip('\n')
		if '>'+alignement_ancestor+'\n' == l:
			a = 1
	
	with open(file_prob,'r') as fin:
		ll = fin.readlines()
	probs = []
	for l in ll:
		if ancestor_file_prob in l:
			l = l.strip('\n').split(',')[2:]
			l = [float(x) for x in l]
			for i,x in enumerate(l):
				if x < th: #arrondir les probabilités sous un certain seuil, considérées impossibles
					l[i] = 0
			l = np.array(l) / sum(l)
			probs.append(l)
	
	dicaa = {0:'A',1:'C',2:'D',3:'E',4:'F',5:'G',6:'H',7:'I',8:'K',9:'L',10:'M',11:'N',12:'P',13:'Q',14:'R',15:'S',16:'T',17:'V',18:'W',19:'Y'}
	
	
	seqs = []
	for k in range(sampling_number):
		seq_random = ''
		for aa,prob in zip(seq,probs):
			if aa != '-' and aa != 'X':
				seq_random += dicaa[choice(range(20),1,p=prob)[0]]
		seqs.append(seq_random)
	
	compteur = {}
	seqsf = []
	
	for s in seqs:
		if s not in seqsf:
			seqsf.append(s)
			compteur[s] = 1
		else:
			compteur[s] +=1
	
	
	with open('ancestor_sampled.fasta','w') as fout:
# 		print(seq)
# 		i = 0
		for s in seqsf:
# 			i+=1
# 			print('>'+str(i))
# 			print(s)
# 			fout.write('>'+str(i)+'\n')
			fout.write('>'+str(compteur[s])+'\n')
			fout.write(s+'\n')
	

#foldx buildmodel séquences alternatives

#conversion en ddG

#histogramme



# score_ref = Repair_model(argv[1],0.1) #utilise la fonction RepairPDB de FoldX sur le modèle initial, jusqu'à ce que l'énergie minimisée soit < th
score_ref = -51.35
#1T9 score_ref = -47.06
# print(score_ref)
# prot_names = list_mut(argv[2],argv[3]) #créer un fichier de mutations faisant le lien antre le modèle initial et les autres séquences mesurées
#1T9prot_names = ['Punidentata1', 'Ppandorae1', 'Pgrasslei2', 'Pmira1', 'Psulfincola2', 'Apompejana2', 'Anc1', 'Anc3', 'Anc4', 'Anc2-T9']
prot_names = ['Punidentata1', 'Ppandorae1', 'Pgrasslei2', 'Pmira1', 'Psulfincola2', 'Apompejana2', 'Anc1-T9', 'Anc3', 'Anc4', 'Anc2-T9']
# print(prot_names)
#Anc1 scores = build_model(argv[1],prot_names,1) #créer le modèle pour les autres protéines mesurées et obtenir le score de stabilité
scores = [-45.05, -40.75, -38.14, -49.18, -38.80, -51.17, -48.54, -52.92, -52.05, -49.86]
#1T9 scores = [-44.31, -40.48, -38.21, -45.89, -35.11, -50.95, -48.83, -53.01, -50.87, -48.91]
# print(scores)

# prot_names.append(argv[2])
# scores.append(score_ref)
# scores = np.array(scores) - score_ref

# a, prot_stabs = optim_corr(scores,argv[4],prot_names) #regression linéaire entre score FoldX et ddG
#MDH Anc1_T9 0.966816915553941 0.23474496907176634 sans fijiensis et unidentata, 55°C
#MDH Anc1    0.8866026290411421 0.26037084347683326 sans fijiensis, 55°C
#SOD Anc1_T9 0.8097394535560207 0.13987495393858615, sans Anc3, Anc6, Ppalmiformis 74°C
#SOD Anc1    0.8707730323001538 0.13585459102041775, sans Anc3, Anc6, Ppalmiformis 74°C
#Hb Anc1 	 0.8569493084947397 0.25744807965908134 sans Anc4, 60°C
#Hb Anc1T9   0.8429714474036163 0.254123859822273 sans Anc4, 60°C

# exit()
########


# sampling_Anc(argv[1],argv[2],argv[3],argv[4],1000,0.01) #nombre de séquences à échantillonner, probabilité pour qu'un résidu soit considéré impossible

# prot = 'Anc1_T9'
# with open(argv[5],'r') as fin: #fichier alignement avec nom séquence ancestrale
# 	ll = fin.readlines()
# sequences = []
# name_seq = []
# w = 1
# for l in ll:
# 	if l[0] != '>' and w == 0:
# 		wt = l.strip('\n')
# 	elif '>'+prot+'\n' == l:
# 		w = 0
# 	elif l[0] == '>' and '>'+prot+'\n' != l:
# 		w = 1
# 
# with open('ancestor_sampled.fasta','a') as fout:
# 	fout.write('>ref\n')
# 	fout.write(wt+'\n')
# 	
# prot_names2 = list_mut('ref', 'ancestor_sampled.fasta')
prot_names2 = [29, 49, 13, 5, 2, 2, 51, 2, 4, 6, 1, 4, 32, 15, 9, 11, 3, 2, 8, 4, 2, 5, 4, 6, 6, 11, 1, 1, 1, 5, 8, 1, 5, 2, 20, 5, 1, 1, 12, 1, 3, 2, 1, 1, 2, 5, 1, 1, 5, 2, 2, 3, 5, 1, 1, 1, 2, 1, 7, 2, 2, 1, 1, 3, 2, 3, 3, 1, 2, 1, 1, 1, 2, 2, 1, 7, 1, 3, 1, 1, 4, 4, 1, 3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
prot_names2 = [str(x) for x in prot_names2]
#1T9 prot_names2 = ['1', '1', '1', '1', '1', '2', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '5', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '5', '1', '1', '1', '1', '1', '1', '1', '4', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '5', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '5', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '4', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '8', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '9', '1', '2', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '5', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '2', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '2', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']
# print(prot_names2)
scores = build_model(argv[6],prot_names2,0)
print(scores)
exit()
scores = [-51.35, -51.43, -50.45, -51.1, -50.51, -50.15, -51.4, -50.53, -51.35, -50.3, -50.44, -50.79, -50.0, -50.89, -50.47, -51.18, -51.38, -49.92, -51.23, -52.61, -51.22, -50.42, -50.81, -51.02, -49.26, -50.24, -49.97, -50.11, -50.67, -50.84, -49.7, -51.38, -51.06, -51.88, -51.9, -51.12, -48.96, -50.78, -51.86, -50.18, -49.14, -50.97, -51.28, -50.51, -48.77, -50.71, -50.99, -49.79, -48.92, -50.32, -49.97, -50.5, -51.84, -51.13, -51.23, -52.32, -49.49, -51.17, -50.88, -50.11, -50.66, -50.5, -49.99, -51.35, -49.93, -50.83, -50.26, -50.76, -51.58, -49.37, -51.28, -50.92, -50.99, -52.75, -51.22, -50.51, -50.31, -49.94]
#1T9scores = [-49.88, -50.63, -49.16, -46.61, -46.23, -46.63, -44.42, -49.7, -48.29, -48.17, -50.82, -49.58, -48.39, -48.42, -46.67, -48.77, -46.89, -49.92, -48.53, -47.14, -48.55, -46.88, -48.2, -47.14, -47.36, -48.81, -48.43, -46.2, -50.75, -49.58, -49.09, -46.21, -49.53, -47.75, -47.71, -48.54, -46.71, -47.63, -47.18, -47.35, -48.91, -49.13, -49.97, -52.33, -46.99, -47.51, -47.29, -47.6, -45.96, -48.08, -46.83, -48.72, -46.62, -47.75, -48.66, -48.0, -47.64, -47.66, -48.17, -47.36, -46.6, -47.55, -47.75, -48.44, -48.77, -50.09, -47.27, -48.67, -47.93, -50.65, -51.25, -47.22, -48.13, -44.21, -48.33, -48.0, -47.49, -48.12, -48.57, -48.3, -48.37, -49.86, -49.08, -47.87, -48.97, -47.05, -47.38, -50.46, -48.6, -47.24, -47.23, -45.41, -47.1, -49.95, -47.7, -50.92, -46.56, -47.89, -46.57, -48.19, -47.88, -48.99, -46.57, -49.03, -46.87, -48.29, -48.46, -50.28, -48.42, -51.55, -51.15, -48.67, -46.89, -47.3, -46.61, -46.87, -49.9, -50.16, -48.6, -47.47, -45.08, -47.09, -48.29, -47.38, -47.9, -47.93, -47.67, -47.85, -47.06, -48.61, -47.77, -48.04, -47.93, -47.63, -48.04, -47.45, -49.05, -48.42, -48.74, -49.7, -47.97, -47.94, -47.53, -47.75, -51.37, -49.02, -47.9, -49.18, -48.04, -48.08, -47.09, -46.34, -48.34, -46.52, -50.51, -48.43, -48.37, -47.13, -48.22, -47.7, -47.98, -46.28, -47.06, -48.81, -48.21, -47.33, -49.6, -47.19, -45.55, -50.64, -50.54, -48.22, -49.23, -47.22, -47.63, -47.33, -45.82, -45.85, -49.64, -47.31, -47.09, -48.68, -48.66, -46.78, -47.92, -47.73, -53.04, -48.19, -47.78, -47.77, -51.6, -49.06, -48.5, -47.65, -48.0, -51.65, -48.26, -48.49, -48.5, -47.91, -46.18, -47.07, -47.98, -48.59, -47.61, -48.29, -48.34, -50.28, -48.03, -46.72, -45.85, -49.05, -48.38, -47.02, -49.7, -49.31, -48.25, -46.79, -48.47, -49.45, -49.75, -49.12, -48.45, -47.72, -48.48, -50.2, -50.65, -48.79, -49.18, -47.41, -48.02, -47.32, -46.33, -48.06, -47.59, -47.54, -48.48, -49.27, -45.59, -47.66, -47.85, -47.95, -49.4, -47.27, -48.32, -49.18, -48.25, -47.24, -50.04, -48.86, -47.98, -48.98, -47.65, -50.14, -48.23, -47.46, -46.49, -46.2, -50.71, -47.11, -48.24, -47.82, -47.86, -47.72, -48.4, -51.03, -48.75, -47.03, -47.91, -51.12, -47.52, -47.99, -44.56, -46.49, -49.11, -48.36, -49.35, -49.18, -47.15, -49.31, -46.74, -47.61, -49.83, -49.44, -44.35, -50.53, -48.66, -47.57, -48.87, -47.01, -49.73, -48.13, -48.24, -48.65, -48.91, -49.86, -47.73, -49.24, -46.68, -49.82, -46.91, -48.06, -49.59, -47.21, -47.1, -48.47, -46.81, -47.86, -48.44, -47.46, -47.56, -48.65, -49.13, -47.59, -47.97, -48.78, -49.8, -47.51, -52.19, -49.3, -46.26, -47.7, -45.48, -49.38, -46.88, -48.84, -47.18, -46.65, -49.45, -48.61, -49.17, -51.46, -47.52, -46.87, -47.72, -48.82, -49.62, -49.35, -47.85, -48.12, -48.83, -47.26, -48.12, -49.19, -47.38, -46.32, -45.07, -47.87, -45.58, -46.16, -47.12, -48.73, -48.96, -50.53, -48.49, -50.48, -50.9, -47.24, -47.49, -47.0, -49.19, -52.07, -46.33, -49.16, -49.13, -43.6, -48.34, -46.89, -49.52, -46.9, -47.15, -46.92, -47.83, -47.69, -48.92, -47.37, -47.59, -49.57, -48.82, -50.36, -45.97, -51.48, -48.19, -51.1, -46.96, -45.93, -48.46, -51.2, -46.56, -49.39, -47.35, -48.39, -46.22, -49.72, -50.12, -47.55, -47.76, -47.78, -49.4, -47.12, -46.25, -47.33, -48.74, -47.09, -48.1, -48.8, -46.83, -46.56, -46.51, -49.04, -45.54, -47.82, -47.76, -47.36, -48.29, -52.27, -47.01, -48.23, -48.51, -47.55, -49.95, -48.19, -46.42, -49.16, -47.5, -48.55, -50.84, -48.14, -49.01, -48.3, -47.32, -48.46, -47.53, -48.44, -46.86, -47.61, -45.69, -46.62, -50.89, -46.93, -48.59, -50.47, -52.06, -47.72, -47.74, -47.08, -48.55, -48.71, -47.17, -48.19, -47.88, -47.72, -48.53, -48.39, -47.68, -47.94, -47.2, -48.43, -47.4, -47.83, -48.54, -49.04, -47.07, -48.28, -48.33, -47.98, -49.66, -47.98, -47.67, -47.18, -49.67, -50.62, -47.19, -49.63, -48.16, -49.27, -47.47, -46.65, -48.32, -48.11, -47.98, -47.79, -47.01, -46.78, -46.94, -48.93, -47.2, -49.44, -47.17, -45.62, -46.61, -46.99, -50.41, -48.21, -51.7, -46.8, -45.24, -47.65, -49.19, -48.07, -48.95, -47.24, -50.07, -49.27, -48.2, -46.98, -47.79, -47.52, -47.4, -48.1, -48.46, -48.26, -47.16, -51.41, -49.33, -49.27, -48.53, -43.96, -46.97, -46.8, -52.48, -48.95, -47.37, -48.57, -47.58, -49.28, -47.85, -46.25, -48.03, -46.34, -48.02, -48.67, -47.5, -47.0, -49.53, -49.31, -46.37, -48.98, -48.16, -47.38, -48.3, -49.4, -47.36, -46.81, -47.94, -48.05, -46.88, -47.04, -48.02, -49.37, -48.52, -46.71, -47.25, -47.29, -50.12, -48.97, -47.18, -46.88, -46.48, -48.82, -47.78, -49.48, -51.12, -50.77, -47.07, -46.91, -46.05, -47.19, -50.5, -47.78, -48.38, -47.49, -47.71, -47.78, -48.92, -46.24, -45.81, -48.4, -46.69, -47.32, -47.66, -50.07, -47.8, -46.53, -48.14, -47.43, -45.9, -47.92, -49.54, -49.6, -49.28, -46.64, -48.3, -48.8, -48.54, -49.48, -48.1, -47.95, -50.24, -47.87, -49.54, -49.29, -47.8, -45.95, -46.97, -46.74, -49.32, -48.12, -47.58, -47.0, -48.63, -45.95, -48.39, -47.3, -48.61, -48.94, -48.02, -49.55, -46.29, -47.15, -48.77, -48.37, -46.5, -47.82, -49.57, -47.58, -52.02, -44.26, -50.61, -48.59, -47.33, -46.56, -48.19, -47.18, -47.78, -48.27, -51.69, -51.18, -49.17, -47.64, -46.39, -45.93, -49.16, -48.26, -47.88, -47.51, -47.1, -48.04, -46.34, -49.93, -48.36, -48.73, -44.98, -47.71, -46.52, -48.1, -45.7, -51.18, -47.78, -46.56, -48.08, -47.44, -46.52, -48.53, -48.41, -47.97, -46.35, -47.82, -47.38, -46.69, -47.36, -47.17, -46.7, -48.13, -47.82, -47.02, -48.91, -47.45, -48.59, -47.08, -49.3, -46.73, -49.83, -50.35, -47.01, -49.21, -47.54, -46.32, -48.56, -47.13, -48.3, -47.95, -47.4, -47.06, -51.87, -48.77, -50.36, -52.0, -46.74, -49.96, -47.88, -46.61, -47.86, -49.68, -47.47, -47.39, -47.47, -48.9, -46.77, -48.18, -48.25, -46.4, -47.19, -49.66, -47.66, -49.47, -47.1, -48.35, -48.18, -47.13, -48.44, -46.63, -47.78, -49.06, -49.38, -48.95, -47.71, -51.87, -46.59, -46.37, -46.55, -51.69, -47.21, -50.76, -49.25, -49.34, -44.36, -50.51, -46.98, -45.99, -50.35, -47.01, -47.64, -48.07, -47.26, -50.87, -46.73, -47.1, -47.8, -46.77, -47.02, -48.02, -49.56, -46.65, -47.19, -47.48, -48.71, -49.97, -48.0, -49.23, -47.1, -47.11, -47.32, -49.64, -47.5, -47.49, -46.88, -47.84, -48.09, -48.45, -48.24, -49.95, -47.93, -47.45, -50.35, -48.29, -46.26, -47.82, -47.98, -47.13, -49.75, -47.64, -50.29, -46.79, -46.14, -49.84, -49.39, -49.19, -49.31, -45.5, -47.68, -47.37, -49.41, -47.44, -47.61, -47.25, -49.27, -48.11, -47.27, -46.26, -49.31, -46.4, -47.22, -48.33, -49.58, -48.23, -46.66, -46.32, -47.61, -49.24, -48.99, -47.42, -48.9, -47.08, -46.3, -47.96, -44.85, -47.26, -49.01, -48.05, -48.0, -46.63, -46.91, -46.98, -48.7, -46.5, -48.21, -46.41, -47.88, -48.23, -49.16, -46.27, -46.61, -47.19, -44.72, -48.79, -47.1, -49.85, -48.14, -50.63, -48.69, -50.54, -47.63, -48.76, -46.58, -49.03, -48.56, -49.81, -47.27, -46.68, -46.85, -50.93, -46.98, -47.64, -47.34, -45.8, -48.26, -46.91, -49.35, -47.57, -46.16, -50.3, -49.49, -48.47, -46.17, -50.49, -48.18, -47.82, -50.02, -47.66, -49.1, -48.23, -46.82, -48.96, -47.71, -48.55, -45.45, -46.39, -49.73, -49.19, -48.27, -50.31, -50.07, -46.27, -47.89, -48.07, -47.51, -47.05, -48.18, -48.36, -48.74, -49.51, -47.41, -48.89, -47.69, -50.82, -47.53, -48.03, -47.28, -47.01, -47.35, -47.48, -45.33, -51.49, -50.7, -48.87, -47.6, -50.35, -46.66]
# print(scores)

scores_weighted = []
for i,k in zip(prot_names2,scores):
	scores_weighted += [(k-score_ref)*a]*int(i)


fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
ax1.hist(scores_weighted,bins=15)
# plt.show()

ax2.scatter(prot_stabs,[0]*len(prot_names))
for p,x in zip(prot_names,prot_stabs):
	ax2.annotate(p,(x,0))
plt.show()



#ligne 170: espèces outlier
#ligne 169: fixer T corrélation
#ligne 74: longueur de la séquence dans modèle 3D
#ligne 169: protéines doublons à ne pas considérer dans la corrélation