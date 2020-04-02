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



niter=1

nd = 3
nran = 500
bs = 1500
msep = (bs**nd/nran)**(1./nd)

if nd == 2:
    rbar = np.sqrt(bs**nd/(np.pi*nran))
if nd == 3:
    rbar = (3.*bs**nd/(4.*np.pi*nran))**(1./3)

#rbar = bs**nd/nran
#rbar = 0.


ranpoints = bs*np.random.uniform(size=(nran,nd))

if nd==2: 
    limits = [[0.,bs],[0.,bs]]#,[0.,bs]]
else:
    limits = [[0.,bs],[0.,bs],[0.,bs]]

if nd==2:
    points1 = ranpoints
else:
    points1 = ranpoints[np.where(ranpoints[:,2]<bs/10)]
    

if nd==2:
    voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])
else:
    voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
    

if nd==3:
    dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
    dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

    print 'Calculating P(k) - Before'    
    real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
    r = FFTPower(real_mesh, mode='1d',dk=0.005)
    Pk1 = r.power

t1 = time.time()

vol1 = [v.get('volume') for v in voro]



for _ in range(niter): 
    #print _,'/',niter-1

    #########
    #RECENTER
    #########    
    meanx=[]
    meany=[]
    if nd==3: meanz=[]
    
    for i in range(nran):
        if nd==2:
            vert=[(a[0],a[1]) for a in voro[i]['vertices']]
        else:
            vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]

        meanx.append( np.mean(np.array(vert)[:,0]) )
        meany.append( np.mean(np.array(vert)[:,1]) )
        if nd==3: meanz.append( np.mean(np.array(vert)[:,2]) )

    if nd==2:
        ranpoints = np.column_stack([meanx,meany])
    else:
        ranpoints = np.column_stack([meanx,meany,meanz])

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

#    if nd==2:
#        voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])
#    else:
#        voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])

    #########
    #VOLUMIZE
    #########
    dr_array =  np.zeros(np.shape(ranpoints))
    for i in range(nran):

        if nd==2:
            dr = rbar - ((voro[i]['volume'])/np.pi)**(1./nd) 
        else:
            dr = rbar - ((3*voro[i]['volume'])/(4*np.pi))**(1./nd)
            
        for face in voro[i]['faces']:
            idx = face['adjacent_cell']
            if idx==i: continue
            dist = ranpoints[idx]-ranpoints[i]
            r_ij = dist/np.linalg.norm(dist)
            dr_array[idx] = dr_array[idx] + dr*r_ij
            #print i,idx,dist,r_ij,dr_array[idx]
 #   mod=0
 #   for i in range(nran):
 #       mod += np.linalg.norm(dr_array[i])
 #   mod /= nran
 #   #print mod
    ranpoints += dr_array

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    if nd==2:
        voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])
    else:
        voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])

t2 = time.time()

print 'Time:',t2-t1

vol2=[v.get('volume') for v in voro]
plt.hist(vol1,bins=20,alpha=.7,label='Original',normed=True)
plt.hist(vol2,bins=20,alpha=.7,label='Recenter+Volumized',normed=True)
plt.xlabel('Voronoi Cell Volume')
plt.legend()
plt.show()

if nd==2:
    points2 = ranpoints

else:
    points2 = ranpoints[np.where(ranpoints[:,2]<bs/10)]

plt.scatter(points2[:,0],points2[:,1],c='r',label='Recenter+Volumized')
plt.scatter(points1[:,0],points1[:,1],alpha=.4,label='Original')
plt.legend()
plt.show()

if nd==3:
    print 'Calculating P(k) - After'   
    dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
    dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

    real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
    r = FFTPower(real_mesh, mode='1d',dk=0.005)
    Pk2 = r.power

    plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
    plt.loglog(Pk2['k'], Pk2['power'].real,label='After')# - Pk.attrs['shotnoise'])
    plt.legend()
    plt.show()

