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
		r[n] = [float(x) for x in l[1:]]
		r[n] = [r[n][0],r[n][2],r[n][3]]
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
	dicnames = {'Punidentata':'Punidentata','Ppandorae':'Ppandorae','Pgrasslei':'Pgrasslei','Ppalmiformis':'Ppalmiformis','Apompejana':'Apompejana','Pfijiensis':'Pfijiensis','Psulfincola':'Psulfincola','Pmira':'Pmira','MDH.T6.asr.0.axes.Anc13':'Anc1-T6','MDH.T6.asr.0.axes.Anc10':'Anc2-T6','MDH.T6.asr.0.axes.Anc7':'Anc3-T6','MDH.T6.asr.0.axes.Anc5':'Anc4-T6','MDH.T6.asr.0.axes.Anc3':'Anc5-T6','MDH.T6.asr.0.axes.Anc6':'Anc6-T6','MDH.T9.asr.0.axes.Anc13':'Anc1-T9'}
	with open(f,'r') as fin:
		ll = fin.readlines()
	r = {}
	for l in ll:
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
	return r



def lin_reg(T,simul,measures,filtre,optim_temp):
	sps = ['Punidentata','Ppandorae','Pgrasslei','Ppalmiformis','Apompejana','Pfijiensis','Psulfincola','Pmira','Anc1-T6','Anc2-T6','Anc3-T6','Anc4-T6','Anc6-T6','Anc1-T9']
	sps2 = ['Pp_diluted' ,'Pp_2' ,'Pf_2' ,'Pf_2_diluted' ,'Anc4_diluted']
	dicsps2 = {'Pp_diluted':'Ppalmiformis' ,'Pp_2':'Ppalmiformis' ,'Pf_2':'Pfijiensis' ,'Pf_2_diluted':'Pfijiensis' ,'Anc4_diluted':'Anc4-T6'}
	
	xdata = [simul[x] for x in sps]
	meas = get_dG(measures,T[0])
	meas = [meas[x] for x in sps+sps2]
# 	ydata = np.array(ydata) - ydata[-1]
	xd = []
	yd = []
	for i,k in enumerate(sps):
		if k not in filtre:
			xd.append(xdata[i])
			yd.append(meas[i])
	for i,k in enumerate(sps2):
		if k not in filtre:
# 			print(k,dicsps2[k],sps.index(dicsps2[k]))
			xd.append(xdata[sps.index(dicsps2[k])])
			yd.append(meas[i+len(sps)])
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
		xdata = [simul[x] for x in sps] + [simul[dicsps2[x]] for x in sps2]
		return r_squared, popt, xdata, meas



def model_calib(measures,predictions,filtre=''):
	
	tempopt = minimize(lin_reg, 25, args=(predictions,measures,filtre,1), method='TNC', bounds=[(10,100)])
	print('Optimal temperature: '+str(tempopt.x[0]))
	
# 	tempopt.x = [50]
	r_square, popt, pred, meas = lin_reg(tempopt.x,predictions,measures,filtre,0)
	
	sps = ['Punidentata','Ppandorae','Pgrasslei','Ppalmiformis','Apompejana','Pfijiensis','Psulfincola','Pmira','Anc1-T6','Anc2-T6','Anc3-T6','Anc4-T6','Anc6-T6','Anc1-T9']
	sps2 = ['Pp_diluted' ,'Pp_2' ,'Pf_2' ,'Pf_2_diluted' ,'Anc4_diluted']
	sps = sps + sps2
	
	plt.scatter(pred,meas,s=10,marker='o')
	for p,x,y in zip(sps,pred,meas):
		plt.annotate(p,(x,y))
	plt.plot([np.min(pred),np.max(pred)],[popt[0]*np.min(pred)+popt[1],popt[0]*np.max(pred)+popt[1]],linewidth = 1)
	print(r_square,popt[0],popt[1])
	plt.show()


measures = get_scores_measures('mesures_MDH.csv')
predictions = get_scores_simul(argv[1])
filtre = ''
# filtre = ['Pfijiensis','Punidentata']
# filtre = 'Pfijiensis'
filtre = ['Pfijiensis','Pf_2','Pf_2_diluted','Ppalmiformis','Pp_2','Anc4-T6']

model_calib(measures,predictions,filtre=filtre)