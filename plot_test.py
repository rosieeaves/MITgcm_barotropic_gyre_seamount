#%%
import numpy as np
import xmitgcm as mit
import matplotlib.pyplot as plt

iter_list = np.linspace(72,21600,300).astype(int).tolist()

data = mit.open_mdsdataset('./run', iters=iter_list)
U = data.U
V = data.V
W = data.W
t = data.time

#%%

print(np.shape(U.mean(dim='time')))
print(np.shape(U.mean(dim='XG').values))
print(np.shape(U.mean(dim='YC').values))

print(np.shape(U.values + V.values))

#%%

KE = np.add(np.square(U.values),np.square(V.values),np.square(W.values))/2
KE_sum = np.sum(np.sum(KE,axis=1),axis=1)
#print(np.shape(KE_sum))

plt.plot(t.values,KE_mean)
plt.show()


# %%
