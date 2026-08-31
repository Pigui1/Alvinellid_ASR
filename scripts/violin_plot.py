import matplotlib.pyplot as plt
import numpy as np

with open('predictions.txt','r') as fin:
	ll = fin.readlines()


def pred_f(x):
	x = -x[1]+x[0]
	coeffs = [0.20773101666932525, -5.003519963824021]
	return coeffs[0]*x + coeffs[1]


def get_scores_measures(f):
	dicnames = {'Pu':'Punidentata','Pk':'Ppandorae','Pg':'Pgrasslei','Pp':'Ppalmiformis','Ap':'Apompejana','Pf':'Pfijiensis','Ps':'Psulfincola','Pm':'Pmira','Anc1':'Anc1-T6','Anc2':'Anc2-T6','Anc3':'Anc3-T6','Anc4':'Anc4-T6','Anc5':'Anc5-T6','Anc6':'Anc6-T6','Anc1-T9':'Anc1-T9','Anc2-T9':'Anc2-T9','Anc3-T9':'Anc3-T9','Anc4-T9':'Anc4-T9','Anc5-T9':'Anc5-T9','Anc6-T9':'Anc6-T9'}
	dicnames2 = {'Pp_diluted':'Pp_diluted' ,'Pp_2':'Pp_2' ,'Pf_2':'Pf_2' ,'Pf_2_diluted':'Pf_2_diluted' ,'Anc4_diluted':'Anc4_diluted'}
	with open(f,'r') as fin:
		ll = fin.readlines()
	r = {}
	for l in ll:
		if l[0] == '#':
			continue
		l = l.strip().split('\t')
		if l[0] in dicnames:
			n = dicnames[l[0]]
		else:
			n = dicnames2[l[0]]
		if -1 not in [float(x) for x in l[1:]]:
			r[n] = [float(x) for x in l[1:]]
			r[n] = [r[n][0],r[n][1],r[n][3]]
			if gibbs_helmotz(0,*r[n]) > 0: #correction for dCp if cold melting above 0°C
				dH = r[n][0]*1000
				dCp = r[n][1]*1000
				Tm = r[n][2] + 273.15
				T = 273.15
				dCp = (dH*(1-T/Tm)) / (Tm - T + T*np.log(T/Tm))
				r[n][1] = dCp/1000
	return r


def gibbs_helmotz(T,dH,dCp,Tm):
	Tm += 273.15
	T += 273.15
	dH *= 1000
	dCp *= 1000
	return -(dH*(1-T/Tm) - dCp*(Tm-T+T*np.log(T/Tm)))/4186


simuls = {'0axis':{'T6':{},'T9':{}},'2axis':{'T6':{},'T9':{}}}

preds = []
for l in ll:
	if l.startswith('>'):
		name = l[1:].strip()
	else:
		preds.append([name,[float(x) for x in l.strip().split(' ')]])

for k,v in preds:
	k = k.split('_')
	for kp in k:
		if '0.axes' in kp:
			dickey = '0axis'
		elif '0missing.model' in kp:
			dickey = '2axis'
		if 'T6' in kp:
			dickey2 = 'T6'
		elif 'T9' in kp:
			dickey2 = 'T9'
		kp = kp.split('.')[-1]
		kp = kp.split('-')
		if kp[0] not in simuls[dickey][dickey2]:
			simuls[dickey][dickey2][kp[0]] = []
		val = pred_f(v)
		for _ in range(int(kp[1])):
			simuls[dickey][dickey2][kp[0]].append(val)

labels = []
maxp = 0

positions_v = {}
positions_model = {'0axis':0,'2axis':0.5}

for c,T in enumerate(['T6','T9']):
	for ax in simuls:
		if ax == '0axis':
			col = 'silver'
		else:
			col='sandybrown'
		for Anc in simuls[ax][T]:
			if Anc not in positions_v:
				positions_v[Anc] = maxp + 2
				if c == 0:
					maxp = positions_v[Anc]
			p = [positions_v[Anc] + positions_model[ax] + c*(maxp + 5)]
			
			parts = plt.violinplot(simuls[ax][T][Anc], positions=p, points=1000, showextrema=False)
			for pc in parts['bodies']:
				pc.set_facecolor(col)
				pc.set_edgecolor('black')
				pc.set_alpha(1)
			
			
			quartile1, medians, quartile3 = np.percentile(simuls[ax][T][Anc], [25, 50, 75])
			whiskers_min, whiskers_max = quartile1, quartile3
			inds = p
			plt.scatter(inds, medians, marker='o', color='white', edgecolors= "mediumvioletred", s=30, zorder=3)
			plt.vlines(inds, quartile1, quartile3, color='mediumvioletred', linestyle='-', lw=1)
			plt.vlines(inds, whiskers_min, whiskers_max, color='mediumvioletred', linestyle='-', lw=1)
			
			plt.annotate(f"{ax} {T} {Anc}",(p[0],-1.5), rotation=90)
	

measures = get_scores_measures('mesures_MDH.csv')
T = 48.504008571138044
x = [0,p[0]]
for k,v in measures.items():
	y = gibbs_helmotz(T,*v)
	if 'Anc' in k:
		if 'T9' in k:
			col = 'dodgerblue'
		else:
			col = 'mediumseagreen'
		plt.scatter(0,[y], c=col, marker='x')
		plt.annotate(k,(0,gibbs_helmotz(T,*v)))
	else:
		plt.plot(x,[y,y], '--', c='k',lw=1)
		plt.annotate(k,(p[0],gibbs_helmotz(T,*v)))

	

plt.show()