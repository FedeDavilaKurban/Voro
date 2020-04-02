import numpy as np 
import recon    
import matplotlib.pyplot as plt
from nbodykit.lab import *
from astropy.table import Table


bs = 100.
nran = 10000
#msep = (bs**3/nran)**(1/3)

np.random.seed(4251241242) 
data=bs*np.random.uniform(size=(nran,3))

dcat = ArrayCatalog(Table(data,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=.07)
Pk1 = r.power



x1=data[:,0]
y1=data[:,1]
z1=data[:,2]

plt.scatter(data[:,0],data[:,1])
plt.show()

rec=recon.Recon(data[:,0],data[:,1],data[:,2],bs)
for i in range(10):
    rec.iterate(i)

rec.apply_shifts_full()
rec.summary()
x,y,z=rec.get_new_radecz(rec.dat)

dcat = ArrayCatalog(Table(np.column_stack((x,y,z)),names=['x','y','z']),BoxSize=bs+20)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=.07)
Pk2 = r.power

plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
plt.loglog(Pk2['k'], Pk2['power'].real,label='After')# - Pk.attrs['shotnoise'])
plt.show()

#print(x[:10]-x1)
x1==x

np.where((x1==x)==False)
np.where((y1==y)==False)
np.where((z1==z)==False)

data[data>bs]-=bs
data[data<0.]+=bs


plt.scatter(x,y)  
plt.show()
