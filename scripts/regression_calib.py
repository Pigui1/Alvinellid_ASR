from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize


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
			r[n] = [r[n][0],r[n][1],r[n][2]]
# 			if gibbs_helmotz(0,*r[n]) > 0: #correction for dCp if cold melting above 0°C
# 				dH = r[n][0]*1000
# 				dCp = r[n][1]*1000
# 				Tm = r[n][2] + 273.15
# 				T = 273.15
# 				dCp = (dH*(1-T/Tm)) / (Tm - T + T*np.log(T/Tm))
# 				r[n][1] = dCp/1000
	return r


def gibbs_helmotz(T,dH,dCp,Tm):
	Tm += 273.15
	T += 273.15
	dH *= 1000
	dCp *= 1000
	return -(dH*(1-T/Tm) - dCp*(Tm-T+T*np.log(T/Tm)))/4186

def get_dG(r,T):
	r_out = {}
	for sp in r:
		r_out[sp] = gibbs_helmotz(T,r[sp][0],r[sp][1],r[sp][2])
	return r_out


def get_scores_simul(f):
	dicnames = {'Punidentata':'Punidentata','Ppandorae':'Ppandorae','Pgrasslei':'Pgrasslei','Ppalmiformis':'Ppalmiformis','Apompejana':'Apompejana','Pfijiensis':'Pfijiensis','Psulfincola':'Psulfincola','Pmira':'Pmira','Anc1':'Anc1-T6','Anc2':'Anc2-T6','Anc3':'Anc3-T6','Anc4':'Anc4-T6','Anc5':'Anc5-T6','Anc6':'Anc6-T6','Anc1_T9':'Anc1-T9','Anc6_T9':'Anc6-T9'}
	with open(f,'r') as fin:
		ll = fin.readlines()
	r = {}
	for l in ll:
		if not l.startswith('#'):
			if l[0] == '>':
				n = dicnames[l[1:].strip()]
			elif n not in r:
				r[n] = [float(x) for x in l.strip().split(' ')]
				#1 of the 2 below, taking WT as reference or not
	# 			r[n] = r[n][0]
				if len(r[n]) > 1:
					r[n] = r[n][0]-r[n][1]
				else:
					r[n] = 0.0
	for _,v in dicnames.items():
		if v not in r:
			r[v] = 0.0
	return r



def lin_reg(T,simul,measures,filtre,optim_temp):
	sps = ['Punidentata','Ppandorae','Pgrasslei','Apompejana','Pfijiensis','Psulfincola','Anc1-T6','Anc3-T6','Anc5-T6','Anc6-T6','Anc1-T9','Anc6-T9']
	
	xdata = [simul[x] for x in sps]
	meas = get_dG(measures,T[0])
	meas = [meas[x] for x in sps]
# 	ydata = np.array(ydata) - ydata[-1]
	xd = []
	yd = []
	for i,k in enumerate(sps):
		if k not in filtre:
			xd.append(xdata[i])
			yd.append(meas[i])
	xd = np.array(xd)
	yd = np.array(yd)
	popt, pcov = curve_fit(lambda x,a,b:a*x+b, xd, yd)
	residuals = yd - popt[0]*xd - popt[1]
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((yd-np.mean(yd))**2)
	r_squared = 1 - (ss_res / ss_tot)
	if optim_temp == 1:
# 		print(1-r_squared)
		return 1-r_squared
	else:
		xdata = [simul[x] for x in sps]
		return r_squared, popt, xdata, meas



def model_calib(measures,predictions,filtre='', display = 1):
	
	tempopt = minimize(lin_reg, 25, args=(predictions,measures,filtre,1), method='TNC', bounds=[(10,100)])
	
# 	tempopt.x = [50]
	r_square, popt, pred, meas = lin_reg(tempopt.x,predictions,measures,filtre,0)
	
	sps = ['Punidentata','Ppandorae','Pgrasslei','Apompejana','Pfijiensis','Psulfincola','Anc1-T6','Anc3-T6','Anc5-T6','Anc6-T6','Anc1-T9','Anc6-T9']
	
	if display == 1:
		plt.scatter(pred,meas,s=10,marker='o')
		for p,x,y in zip(sps,pred,meas):
			plt.annotate(p,(x,y))
		plt.plot([np.min(pred),np.max(pred)],[popt[0]*np.min(pred)+popt[1],popt[0]*np.max(pred)+popt[1]],linewidth = 1,label=str(r_square))
		plt.legend()
		print('Optimal temperature: '+str(tempopt.x[0]))
		print(r_square,popt[0],popt[1])
		plt.savefig('figure.svg')
		plt.show()
	
	return r_square


measures = get_scores_measures('mesures_SOD.csv')
predictions = get_scores_simul(argv[1])
filtre = ''
# filtre = ['Pfijiensis','Punidentata']
# filtre = 'Pfijiensis'
# filtre = ['Pfijiensis','Pf_2','Pf_2_diluted','Ppalmiformis','Pp_2','Anc4-T6']
# filtre = ['Pfijiensis','Pf_2','Ppalmiformis','Pp_2']
filtre = []

r_square = model_calib(measures,predictions,filtre=filtre, display=0)
print(f"initial r_square {r_square}")

rmv_max = 2
r_square_opt = 0
spc_removed = None
for i in range(rmv_max):
	for k in measures:
		if k not in filtre:
			filtre2 = filtre + [k]
			r_square = model_calib(measures,predictions,filtre=filtre2, display=0)
			if r_square > r_square_opt:
				spc_removed = k
				r_square_opt = r_square
	filtre = filtre + [spc_removed]
	print(f"{spc_removed} removed, r_square {r_square_opt}")

# _ = model_calib(measures,predictions,filtre=filtre, display=1)


