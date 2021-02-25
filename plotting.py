#%%
# load data
import numpy as np
import xmitgcm as mit
import matplotlib.pyplot as plt

class Data:

    def __init__(self,dataDir,iter_list,geom,dt,dumpFreq,f0,beta):
        self.dataDir = dataDir
        self.iter_list = iter_list
        self.data = mit.open_mdsdataset(data_dir=dataDir, iters=iter_list, geometry=geom)
        
        self.dt = dt
        self.dumpFreq = dumpFreq
        self.f0 = f0
        self.beta = beta

        self.Nx = len(self.data.XC)
        self.Ny = len(self.data.YC)
        self.dx = self.data.dxG.values[0][0] # works if the spacing is constant
        self.dy = self.data.dyG.values[0][0] # works if the spacing is constant
        self.x = self.data.XC.values
        self.y = self.data.YC.values

        self.volume = np.mean(self.data.hFacC * self.data.drF * self.data.rA,2)
        self.t = self.data.time
        self.time = ((self.t.values/(10**9))*self.dt)
        self.f = [np.ones(self.Nx)*(self.f0+(i+0.5)*self.dy*self.beta) for i in range(self.Ny)]
        
        self.PH = self.data.PH
        self.Depth = self.data.Depth
        self.U = self.data.U
        self.V = self.data.V
        self.W = self.data.W
        self.T = self.data.T

    def plot_streamU(self,times):
        
        h = self.Depth.values 
        U = self.U.values 

        for t in times:
            time_index = ((t*86400)/self.dumpFreq - 1).astype(int)
            streamU = self.dy*(np.multiply(h,U[time_index]))

            X,Y = np.meshgrid(self.x,self.y)
            plt.contour(X,Y,streamU)
            plt.colorbar()
            plt.show()

    def plot_streamV(self,times):
        
        h = self.Depth.values 
        V = self.V.values 

        for t in times:
            time_index = ((t*86400)/self.dumpFreq - 1).astype(int)
            streamV = self.dy*(np.multiply(h,V[time_index]))

            X,Y = np.meshgrid(self.x,self.y)
            plt.contour(X,Y,streamV)
            plt.colorbar()
            plt.show()

    def plot_bathymetry(self):
        
        X,Y = np.meshgrid(self.x,self.y)
        plt.contourf(X,Y,self.Depth.values)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        cbar = plt.colorbar()
        cbar.set_label('Depth (m)')
        plt.show()

    def plot_KE_volAvg(self):

        KE = (np.add(np.square(self.U.values),np.square(self.V.values),np.square(self.W.values))/2)*self.volume.values
        KE_sum = np.sum(np.sum(KE,axis=1),axis=1)
        KE_volAvg = KE_sum/np.sum(self.volume.values)

        plt.plot(self.time,KE_volAvg)
        plt.xlabel('time (days)')
        plt.ylabel('Volume averaged KE')
        plt.show()

    def plot_f(self):
        
        X,Y = np.meshgrid(self.x,self.y)
        plt.contourf(X,Y,self.f)
        plt.colorbar()
        plt.show()

    def plot_vorticity(self,times):

        u = self.U.values
        v = self.V.values

        for t in times:

            time_index = ((t*86400)/self.dumpFreq - 1).astype(int)

            relVort = np.zeros((self.Nx,self.Ny))

            # calculate vorticity matrix

            for i in range(self.Nx):
                for j in range(self.Ny):

                    #dvdx
                    if i==0:
                        dvdx = (v[time_index][i+1][j] - v[time_index][i][j])/self.dx
                    elif i==self.Nx-1:
                        dvdx = (v[time_index][i][j] - v[time_index][i-1][j])/self.dx
                    else:
                        dvdx = (v[time_index][i+1][j] - v[time_index][i-1][j])/(2*self.dx)

                    #dudy
                    if j==0:
                        dudy = (u[time_index][i][j+1] - u[time_index][i][j])/self.dy
                    elif j==self.Ny-1:
                        dudy = (u[time_index][i][j] - u[time_index][i][j-1])/self.dy
                    else:
                        dudy = (u[time_index][i][j+1] - u[time_index][i][j-1])/(2*self.dy)

                    relVort[i][j] = dvdx - dudy

            vorticity = relVort

            # plot

            X,Y = np.meshgrid(self.x,self.y)
            plt.contour(X,Y,vorticity)
            plt.xlabel('x (km)')
            plt.ylabel('y (km)')
            cbar = plt.colorbar()
            cbar.set_label('Vorticity (s$^{-1}$)')
            plt.title('Day ' + str((self.time[time_index]/86400).astype(int)))
            plt.show()

    def plot_T(self,times):

        x_noBoundaries = self.x[1:self.Nx-1]
        y_noBoundaries = self.y[1:self.Ny-1]

        T_noBoundaries = np.delete(np.delete(np.delete(np.delete(self.T.values,self.Nx-1,1),0,1),self.Ny-1,2),0,2)

        X_noBoundaries,Y_noBoundaries = np.meshgrid(x_noBoundaries,y_noBoundaries)

        for t in times:
            plt.contourf(X_noBoundaries,Y_noBoundaries,T_noBoundaries[t])
            plt.colorbar()
            plt.show()

        


#%%       
#iters = np.linspace(72,21600,300).astype(int).tolist()
iters = np.linspace(160,58400,365).astype(int)
baro_wind_mount = Data(dataDir='./run',iter_list=iters,geom='cartesian',dt=5400,dumpFreq=864000,f0=0.7*10**(-4),beta=2*10**(-11))

#%%
t = np.linspace(300,3000,10).astype(int) 
baro_wind_mount.plot_streamU(times=t)

#%%

print(baro_wind_mount.x)

# %%

#t = np.linspace(30,300,10).astype(int) 
baro_wind_mount.plot_vorticity(times=t)

# %%
