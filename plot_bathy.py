#%%

import numpy as np
import matplotlib.pyplot as plt 

bathy = np.load('bathy.npy')
print(bathy)

x = np.linspace(1,62,62)*20
y = np.linspace(1,62,62)*20
X,Y = np.meshgrid(x,y)

plt.contourf(X,Y,bathy)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
cbar = plt.colorbar()
cbar.set_label('Depth (m)')
plt.show()

# %%
