#%%

import numpy as np
import xmitgcm as mit
import matplotlib.pyplot as plt 

Nx = 62
Ny = 62

x = np.linspace(1,Nx,Nx)
y = np.linspace(1,Ny,Ny)

L = Nx
W = Ny

H = 5000 # depth of basin
h = 500 # height of seamount

bathy = np.zeros((Nx,Ny))

for i in range(Nx):
    for j in range(Ny):
        if i!=0 and i!=Nx-1 and j!=0 and j!=Ny-1:
            bathy[i][j] = -(H-h*np.sin((np.pi*x[i])/L)*np.sin((np.pi*y[j])/W)) 


X,Y = np.meshgrid(x,y)

#%%

plt.contourf(X,Y,bathy)
plt.colorbar()
plt.show()

np.save('bathy', bathy)

# %%

mit.utils.write_to_binary(bathy.flatten(), 'bathy.bin')

# %%
