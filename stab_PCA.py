import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv('dG_site_feature_Fig3.csv', names=['ddg_Q','ddg_E','ddg_N','ddg_H','ddg_D','ddg_R','ddg_K','ddg_T','ddg_S','ddg_A','ddg_G','ddg_M','ddg_C','ddg_L','ddg_V','ddg_I','ddg_W','ddg_Y','ddg_F','ddg_P'], header=0, sep=";")


features = ['ddg_Q','ddg_E','ddg_N','ddg_H','ddg_D','ddg_R','ddg_K','ddg_T','ddg_S','ddg_A','ddg_G','ddg_M','ddg_C','ddg_L','ddg_V','ddg_I','ddg_W','ddg_Y','ddg_F','ddg_P']

amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# Separating out the features
x = df.loc[:, features].values
for i,k in enumerate(x):
	x[i] = k - np.sum(k)/20 #mean
# 	x[i] = k - k[7] #arbitrary aa reference
# 	x[i] = k - np.median(k) #median
# 	x[i] = k - np.max(k) #most stable aa


# Standardizing the features
#x = StandardScaler().fit_transform(x)


pca = PCA()

principalComponents = pca.fit_transform(x)

col_names = ['principal component '+str(k) for k in range(1,21)]
principalDf = pd.DataFrame(data = principalComponents, columns = col_names)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)


# for k in principalDf.loc:
# 	ax.scatter([0,0], [k[0],k[1]])
# ax.grid()

for k in zip(pca.components_[0],pca.components_[1]):
	ax.plot([0,k[0]], [0,k[1]])
ax.grid()

#print(pca.components_)
print(pca.explained_variance_/np.sum(pca.explained_variance_))

# plt.show()

# axs = range(20)
# print(features)
# for ax in axs:
# 	print(pca.components_[ax])
# 	axis = pca.components_[ax] - np.sum(pca.components_[ax])/20
# 	sorted_axis = axis.argsort()
# 	print(sorted_axis)
# 	features2 = [features[x] for x in sorted_axis]
# 	plt.plot(range(20),sorted(axis))
# 	print(features2)


#plt.show()


dic_features = {x:features.index('ddg_'+x) for x in amino_acids}

with open('thermostable.txt','w') as fout:
	fout.write('\t'.join(amino_acids)+'\n')
	for v in pca.components_:
		v = [str(v[dic_features[x]]) for x in amino_acids]
		fout.write('\t'.join(v)+'\n')




