#%%
# generate tau
import xmitgcm as mit
import numpy as np 
import matplotlib.pyplot as plt 

Nx=248
Ny=248 
dx = 5
dy = 5
L = Nx*dx
W = Ny*dy

tauMax=0.1
x=[(i+0.5)*dx for i in range(Nx)]
y=[(i+0.5)*dy for i in range(Ny)]

tau = np.zeros((Nx,Ny))
for i in range(Nx):
    for j in range(Ny):
        tau[i][j] = -tauMax*np.cos((np.pi*y[j])/W)

tau = np.transpose(tau)

#%%
# plot tau

X,Y = np.meshgrid(x,y)
plt.contourf(X,Y,tau)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
cbar = plt.colorbar()
cbar.set_label = 'Wind stress (Nm$^{-2}$)'
plt.show()


# %%
# generate files 

np.save('wind',tau)

mit.utils.write_to_binary(tau.flatten(), 'wind.bin')

# %%
