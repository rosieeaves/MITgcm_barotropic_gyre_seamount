#%%
import numpy as np
import xmitgcm as mit
import matplotlib.pyplot as plt

iter_list = np.linspace(72,21600,300).astype(int).tolist()

data = mit.open_mdsdataset('./run', iters=iter_list, geometry='cartesian')
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

print(np.shape(data.Depth.values))

x = np.linspace(1,62,62)*20
y = np.linspace(1,62,62)*20
X,Y = np.meshgrid(x,y)

plt.contourf(X,Y,data.Depth.values)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
cbar = plt.colorbar()
cbar.set_label('Depth (m)')
plt.show()

#%%

volume = np.mean(data.hFacC * data.drF * data.rA,2)

print(np.shape(volume.values))

plt.contourf(X,Y,volume.values)
plt.colorbar()
plt.show()

#%%

KE = np.add(np.square(U.values),np.square(V.values),np.square(W.values))/2
KE_sum = np.sum(np.sum(KE,axis=1),axis=1)
#print(np.shape(KE_sum))

plt.plot(t.values,KE_mean)
plt.show()


# %%

