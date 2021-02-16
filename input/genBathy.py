#%%
# generate bathy

import numpy as np
import xmitgcm as mit
import matplotlib.pyplot as plt 

Nx = 248
Ny = 248
dx = 5
dy = 5

x=[(i+0.5)*dx for i in range(Nx)]
y=[(i+0.5)*dy for i in range(Ny)]

L = Nx*dx
W = Ny*dy

H = 5000 # depth of basin
h = 500 # height of seamount

bathy = np.zeros((Nx,Ny))

for i in range(Nx):
    for j in range(Ny):
        if i!=0 and i!=Nx-1 and j!=0 and j!=Ny-1:
            bathy[i][j] = -(H-h*np.sin((np.pi*x[i])/L)*np.sin((np.pi*y[j])/W)) 

bathy = np.transpose(bathy)

#%%
# plot bathy

X,Y = np.meshgrid(x,y)
plt.contourf(X,Y,bathy)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
cbar = plt.colorbar()
cbar.set_label('Depth (m)')
plt.show()

# %%
# create bathy.bin

np.save('bathy', bathy)

mit.utils.write_to_binary(bathy.flatten(), 'bathy.bin')

# %%
