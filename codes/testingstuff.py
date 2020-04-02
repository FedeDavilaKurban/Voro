"""
Attempting to do voro.py with a way of equalizing the volumes
+
Recenter

v2: 2d and 3d in one code
"""
import numpy as np
import matplotlib.pyplot as plt
import pyvoro
from nbodykit.lab import *
import ascii
from ascii import Table


niter=25

nd = 3
nran = 16**3
bs = 1500
msep = (bs**nd/nran)**(1./nd)
smooth = 1.3

rbar = (3.*bs**nd/(4.*np.pi*nran))**(1./3)

#rbar = bs**nd/nran
#rbar = 0.

np.random.seed(4294967287)
ranpoints = bs*np.random.uniform(size=(nran,nd))

limits = [[0.,bs],[0.,bs],[0.,bs]]

points1 = ranpoints[np.where(ranpoints[:,2]<bs/10)]
    

voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
    

dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - Before'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.005)
Pk1 = r.power

#t1 = time.time()

vol1 = [v.get('volume') for v in voro]
norm = np.zeros(niter)
for _ in range(niter): 
    print _,'/',niter-1

    #########
    #RECENTER
    #########    
    meanx=[]
    meany=[]
    meanz=[]
    
    for i in range(nran):
        vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]

        meanx.append( np.mean(np.array(vert)[:,0]) )
        meany.append( np.mean(np.array(vert)[:,1]) )
        meanz.append( np.mean(np.array(vert)[:,2]) )

    ranpoints = np.column_stack([meanx,meany,meanz])

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    #########
    #VOLUMIZE
    #########
    dr_array =  np.zeros(np.shape(ranpoints))
    for i in range(nran):

        dr = rbar - ((3*voro[i]['volume'])/(4*np.pi))**(1./nd)
            
        for face in voro[i]['faces']:
            idx = face['adjacent_cell']
        #for testing in range(20):
        #    idx = np.random.randint(0,nran)
            if idx==i: continue
            dist = ranpoints[idx]-ranpoints[i]
            r_ij = dist/np.linalg.norm(dist)
            dr_array[idx] = dr_array[idx] + smooth*dr*r_ij
            #print i,idx,dist,r_ij,dr_array[idx]
    
    norm[_] = np.mean(np.linalg.norm(dr_array,axis=1))
    
    ranpoints += dr_array

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])

#    if _ in [10,20,30,40,50,60,70,80,90]:
#        dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
#        dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

#        real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
#        r = FFTPower(real_mesh, mode='1d',dk=0.005)
#        Pk = r.power

#        plt.loglog(Pk['k'], Pk['power'].real,label='niter={}'.format(_))# - Pk.attrs['shotnoise'])


#t2 = time.time()

#ascii.write(ranpoints,'volumized_cat/volumized_{}_{}.txt'.format(nran,niter),format='no_header',overwrite=True)

#print 'Time:',t2-t1

vol2=[v.get('volume') for v in voro]

################
###
#Calculate Voronoi Volumes in CCVT for Comparison (MPE Cluster)
###

rcat = CSVCatalog('../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
rcat['Position']*=bs
real_mesh = rcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[bs,bs,bs],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.005)
Pkccvt = r.power

ccvt_points = rcat['Position'].compute()
ccvt_voro = pyvoro.compute_voronoi(ccvt_points,limits,msep,periodic=[True,True,True])
ccvt_vol=[v.get('volume') for v in ccvt_voro]

################
###
#Calculate Voronoi Volumes in Glass for Comparison (MPE Cluster)
###

gcat = CSVCatalog('../data/glass_16.dat',names=['x','y','z'])
gcat['Position'] = transform.StackColumns(gcat['x'], gcat['y'], gcat['z'])
gcat['Position']*=bs
real_mesh = gcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[bs,bs,bs],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.005)
Pkglass = r.power

gpoints = gcat['Position'].compute()
gvoro = pyvoro.compute_voronoi(gpoints,limits,1.5*msep,periodic=[True,True,True])
gvol=[v.get('volume') for v in gvoro]


plt.hist(vol1,bins=20,alpha=.7,label='Original',normed=True)
plt.hist(vol2,bins=20,alpha=.7,label='Recenter+Volumized',normed=True)
plt.hist(gvol,bins=20,alpha=.7,label='Glass',normed=True)
plt.hist(ccvt_vol,bins=20,alpha=.7,label='CCVT',normed=True)
plt.xlabel('Voronoi Cell Volume')
plt.legend()
plt.show()

points2 = ranpoints[np.where(ranpoints[:,2]<bs/10)]

plt.scatter(points2[:,0],points2[:,1],c='r',label='Recenter+Volumized')
plt.scatter(points1[:,0],points1[:,1],alpha=.4,label='Original')
plt.legend()
plt.show()

print 'Calculating P(k) - After'   
dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.005)
Pk2 = r.power

plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
plt.loglog(Pk2['k'], Pk2['power'].real,label='After')# - Pk.attrs['shotnoise'])
plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
plt.loglog(Pk2['k'],((bs**3/nran)*(Pk2['k']*msep/(2.*np.pi))**4),color='k',ls=':') 
plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass')

plt.hlines((bs**nd)/nran, Pk1['k'][0],Pk2['k'][-1], colors='k',linestyles=':')
plt.vlines(2.*np.pi/msep,10,10E9,colors='k',linestyles=':')
plt.legend()
plt.show()

plt.plot(range(niter),norm,label=r'$\Delta r$')
plt.yscale('log')
#plt.xscale('log')
plt.legend()
plt.show()
