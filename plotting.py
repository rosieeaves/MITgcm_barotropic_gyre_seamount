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
        self.x = np.linspace(1,self.Nx,self.Nx)*self.dx/1000
        self.y = np.linspace(1,self.Ny,self.Ny)*self.dy/1000

        self.volume = np.mean(self.data.hFacC * self.data.drF * self.data.rA,2)
        self.t = self.data.time
        self.time = ((self.t.values/(10**9))*self.dt)/self.dumpFreq
        self.f = [np.ones(self.Nx)*(self.f0+(i+0.5)*self.dy*self.beta) for i in range(self.Ny)]
        
        self.PH = self.data.PH
        self.Depth = self.data.Depth
        self.U = self.data.U
        self.V = self.data.V
        self.W = self.data.W

    def plot_streamfunction(self,times):
        x_noBoundaries = self.x[1:self.Nx-1]
        y_noBoundaries = self.y[1:self.Ny-1]
        f_noBoundaries = np.delete(np.delete(np.delete(np.delete(self.f,self.Nx-1,0),0,0),self.Ny-1,1),0,1)

        PH_noBoundaries = np.delete(np.delete(np.delete(np.delete(self.PH.values,self.Nx-1,1),0,1),self.Ny-1,2),0,2)

        X_noBoundaries,Y_noBoundaries = np.meshgrid(x_noBoundaries,y_noBoundaries)

        for t in times:
            streamfunc = np.divide(PH_noBoundaries[t-1],f_noBoundaries) 

            plt.contour(X_noBoundaries,Y_noBoundaries,streamfunc)
            plt.xlabel('x (km)')
            plt.ylabel('y (km)')
            cbar = plt.colorbar()
            cbar.set_label('Streamfunction')
            plt.title('Day' + str(t-1))
            plt.show()

    def plot_bathymetry(self):
        
        X,Y = np.meshgrid(self.x,self.y)
        plt.contourf(X,Y,self.Depth.values)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        cbar = plt.colorbar()
        cbar.set_label('Depth (m)')
        plt.show()

    def plot_volume(self):

        X,Y = np.meshgrid(self.x,self.y)
        plt.contourf(X,Y,self.volume.values)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        cbar = plt.colorbar()
        cbar.set_label('Volume (m$^2$)')
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


#%%       
iters = np.linspace(72,21600,300).astype(int).tolist()
baro_wind_mount = Data(dataDir='./run',iter_list=iters,geom='cartesian',dt=1200,dumpFreq=86400,f0=10**(-4),beta=10**(-11))

#%%

t = np.linspace(6,31,6).astype(int)
baro_wind_mount.plot_streamfunction(times=t)


# %%

baro_wind_mount.plot_KE_volAvg()

# %%
