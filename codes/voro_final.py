"""
Attempting to do voro.py with a way of equalizing the volumes
+
Recenter
"""



###
### 2D
###

import numpy as np
import matplotlib.pyplot as plt
import pyvoro
niter=2points00

nd = 2
nran = 100
bs = 1.
msep = (bs**nd/nran)**(1./nd)

ranpoints = bs*np.random.uniform(size=(nran,nd))
limits = [[0.,bs],[0.,bs]]#,[0.,bs]]

points1 = ranpoints
plt.scatter(points1[:,0],points1[:,1])
plt.show()

voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])#,True])

#################
#Calculo del Pk
################
#dcat = ArrayCatalog(Table(ranpoints,names=['x','y']),BoxSize=bs)
#dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'])

#print 'Calculating P(k) - Before'    
#real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
#r = FFTPower(real_mesh, mode='1d',dk=0.05)
#Pk1 = r.power


vol=[v.get('volume') for v in voro]
#print len(vol)
plt.hist(vol,bins=20,alpha=.7,label='Original',normed=True)
plt.xlabel('Voronoi Cell Volume')
plt.legend()

#mshift=[]

for _ in range(niter): 
    print _
    #voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])#,True])
    #########
    #RECENTER
    #########    
#    meanx=[]
#    meany=[]
    #meanz=[]
    #shift=[]
#    for i in range(nran):
#        vert=[(a[0],a[1]) for a in voro[i]['vertices']]
#        meanx.append( np.mean(np.array(vert)[:,0]) )
#        meany.append( np.mean(np.array(vert)[:,1]) )
        #meanz.append( np.mean(np.array(vert)[:,2]) )
        #x,y = voro[i]['original']
        #shift.append( np.sqrt((x-meanx[-1])**2+(y-meany[-1])**2) ) 
    #mshift.append( np.mean(np.array(shift)) )
    
#    ranpoints = np.column_stack([meanx,meany])


    #########
    #VOLUMIZE
    #########
    rbar=((bs**2)/(np.pi*nran))**(1./2)
    #rbar = msep
    #print 'A'
    dr_array =  np.zeros(np.shape(ranpoints))
    for i in range(nran):
        #print 'B'
        dr = rbar - ((voro[i]['volume'])/(np.pi*nran))**(1./2)

        for face in voro[i]['faces']:
            idx = face['adjacent_cell']
            if idx==i: continue
            dist = ranpoints[i]-ranpoints[idx]
            r_ij = dist/np.linalg.norm(dist)
            dr_array[idx] = dr_array[idx] + dr*r_ij
            #print i,idx,dist,r_ij,dr_array[idx]
        #print 'C'
    ranpoints += dr_array

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])

#voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])
vol=[v.get('volume') for v in voro]
#print len(vol)
plt.hist(vol,bins=20,alpha=.7,label='Volumized',normed=True)
plt.xlabel('Voronoi Cell Volume')
plt.legend()
plt.show()

points2 = ranpoints
plt.scatter(points1[:,0],points1[:,1])
plt.scatter(points2[:,0],points2[:,1],alpha=.7,c='r')
plt.show()

######################################################################
######################################################################
######################################################################
###
### 3D
###



#################################
#Equalize Volume v1 + Recenter
#################################
import numpy as np
import matplotlib.pyplot as plt
import pyvoro
niter=100

nd = 3
nran = 100
bs = 1.
msep = (bs**nd/nran)**(1./nd)

ranpoints = bs*np.random.uniform(size=(nran,nd))
limits = [[0.,bs],[0.,bs],[0.,bs]]

points1 = ranpoints[np.where(ranpoints[:,2]<.1)]

voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])

#################
#Calculo del Pk
################
dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - Before'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.05)
Pk1 = r.power


vol=[v.get('volume') for v in voro]
plt.hist(vol,bins=20,alpha=.7,label='Original')
plt.xlabel('Voronoi Cell Volume')
plt.legend()

#mshift=[]

for _ in range(niter): 
    print _

    #########
    #RECENTER
    #########    
    meanx=[]
    meany=[]
    meanz=[]
    shift=[]
    for i in range(nran):
        vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]
        meanx.append( np.mean(np.array(vert)[:,0]) )
        meany.append( np.mean(np.array(vert)[:,1]) )
        meanz.append( np.mean(np.array(vert)[:,2]) )
        x,y,z = voro[i]['original']
        shift.append( np.sqrt((x-meanx[-1])**2+(y-meany[-1])**2+(z-meanz[-1])**2) ) 
    #mshift.append( np.mean(np.array(shift)) )
    
    ranpoints = np.column_stack([meanx,meany,meanz])


    #########
    #VOLUMIZE
    #########
    rbar=((3*bs**3)/(4*np.pi*nran))**(1./3)
    dr_array =  np.zeros(np.shape(ranpoints))

    for i in range(nran):
        #print i
        dr = rbar - ((3*voro[i]['volume'])/(4*np.pi*nran))**(1./3)

        for face in voro[i]['faces']:
            idx = face['adjacent_cell']
            dist = ranpoints[i]-ranpoints[idx]
            r_ij = dist/np.linalg.norm(dist)
            dr_array[idx] = dr_array[idx] + dr*r_ij
        
    ranpoints += dr_array

    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs


voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
vol=[v.get('volume') for v in voro]
plt.hist(vol,bins=20,alpha=.7,label='Volumized')
plt.xlabel('Voronoi Cell Volume')
plt.legend()
plt.show()

points2 = ranpoints[np.where(ranpoints[:,2]<.1)]
plt.scatter(points1[:,0],points1[:,1])
plt.scatter(points2[:,0],points2[:,1],c='r')
plt.show()

dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - Before'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.05)
Pk2 = r.power
plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
plt.loglog(Pk2['k'], Pk2['power'].real,label='After')# - Pk.attrs['shotnoise'])
plt.legend()
plt.show()



