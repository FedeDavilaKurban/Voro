"""
Attempting to do voro.py with a way of equalizing the volumes
+
Recenter

v2: 2d and 3d in one code
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import pyvoro
from nbodykit.lab import *
import time

############
def copybox(pos,buf,box):

    g1 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box and box-buf<i[2]<box])
    g1[:, 0] -= box
    g1[:, 1] -= box
    g1[:, 2] -= box
    g2 = np.asarray([i for i in pos if box-buf<i[1]<box and box-buf<i[2]<box])
    g2[:, 1] -= box
    g2[:, 2] -= box
    g3 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box and box-buf<i[2]<box])
    g3[:, 0] += box
    g3[:, 1] -= box
    g3[:, 2] -= box
    g4 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[2]<box])
    g4[:, 0] -= box
    g4[:, 2] -= box
    g5 = np.asarray([i for i in pos if box-buf<i[2]<box])
    g5[:, 2] -= box
    g6 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[2]<box])
    g6[:, 0] += box
    g6[:, 2] -= box
    g7 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf and box-buf<i[2]<box])
    g7[:, 0] -= box
    g7[:, 1] += box
    g7[:, 2] -= box

    g8 = np.asarray([i for i in pos if 0<i[1]<buf and box-buf<i[2]<box])
    g8[:, 1] += box
    g8[:, 2] -= box
    g9 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf and box-buf<i[2]<box])
    g9[:, 0] += box
    g9[:, 1] += box
    g9[:, 2] -= box
    g10 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box])
    g10[:, 0] -= box
    g10[:, 1] -= box
    g11 = np.asarray([i for i in pos if box-buf<i[1]<box])
    g11[:, 1] -= box
    g12 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box])
    g12[:, 0] += box
    g12[:, 1] -= box
    g13 = np.asarray([i for i in pos if box-buf<i[0]<box])
    g13[:, 0] -= box
    g14 = np.asarray([i for i in pos if 0<i[0]<buf])
    g14[:, 0] += box
    g15 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf])
    g15[:, 0] -= box
    g15[:, 1] += box
    g16 = np.asarray([i for i in pos if 0<i[1]<buf])
    g16[:, 1] += box
    g17 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf])
    g17[:, 0] += box
    g17[:, 1] += box
    g18 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box and 0<i[2]<buf])
    g18[:, 0] -= box
    g18[:, 1] -= box
    g18[:, 2] += box
    g19 = np.asarray([i for i in pos if box-buf<i[1]<box and 0<i[2]<buf])
    g19[:, 1] -= box
    g19[:, 2] += box
    g20 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box and 0<i[2]<buf])
    g20[:, 0] += box
    g20[:, 1] -= box
    g20[:, 2] += box
    g21 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[2]<buf])
    g21[:, 0] -= box
    g21[:, 2] += box
    g22 = np.asarray([i for i in pos if 0<i[2]<buf])
    g22[:, 2] += box
    g23 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[2]<buf])
    g23[:, 0] += box
    g23[:, 2] += box
    g24 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf and 0<i[2]<buf])
    g24[:, 0] -= box
    g24[:, 1] += box
    g24[:, 2] += box
    g25 = np.asarray([i for i in pos if 0<i[1]<buf and 0<i[2]<buf])
    g25[:, 1] += box
    g25[:, 2] += box
    g26 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf and 0<i[2]<buf])
    g26[:, 0] += box
    g26[:, 1] += box
    g26[:, 2] += box

    pos = np.vstack([pos, g1, g2, g3, g4, g5, g6, g7, g8, g9,
                    g10, g11, g12, g13, g14, g15, g16, g17,
                    g18, g19, g20, g21, g22, g23, g24, g25, g26])
    return pos
###########


recenter = True
volumize = True
periodicVoro = False
smooth = .1 #se usa en el volumizado 
niter=10
buf=150.

nd = 3
nran = 16**3
bs = 1500.
msep_real = (bs**nd/nran)**(1./nd)

dk = 0.007

#TEST multiplico la separacion media para usarla en la teselacion Voronoi
#msep = np.sqrt(2)*msep_real
msep = 1.2*msep_real
#msep = 2.*msep_real

rbar = (3.*bs**nd/(4.*np.pi*nran))**(1./3)

#rbar = bs**nd/nran
#rbar = 0.
#rbar = msep

np.random.seed(2022342)
ranpoints = bs*np.random.uniform(size=(nran,nd))

##USO EL RBAR EN VEZ DEL MSEP!!!
rbar = (3.*bs**nd/(4.*np.pi*len(ranpoints)))**(1./3)
msep = 1.2*rbar

####Calculate Pk####################
print 'Calculating P(k) - Before'    
dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pk1 = r.power

vol1 = [v.get('volume') for v in voro]
####################################

limits = [[0.,bs],[0.,bs],[0.,bs]]
limits_cp = [[-150,bs+buf],[-150,bs+buf],[-150,bs+buf]]

points1 = ranpoints[np.where(ranpoints[:,2]<bs/10)]

ranpoints = copybox(ranpoints,buf,bs)
msep_real = ((bs+buf)**nd/len(ranpoints))**(1./nd)
msep = 1.2*msep_real

##USO EL RBAR EN VEZ DEL MSEP!!!
rbar = (3.*bs**nd/(4.*np.pi*len(ranpoints)))**(1./3)
msep = 1.2*rbar


voro = pyvoro.compute_voronoi(ranpoints,limits_cp,msep)



#test
#neg=0
#pos=0

for _ in range(niter): 
    print _,'/',niter-1
    if recenter==True:
        #########
        #RECENTER
        #########    
        meanx=[]
        meany=[]
        meanz=[]

        for i in range(len(ranpoints)):
            vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]

            meanx.append( np.mean(np.array(vert)[:,0]) )
            meany.append( np.mean(np.array(vert)[:,1]) )
            meanz.append( np.mean(np.array(vert)[:,2]) )

        ranpoints = np.column_stack([meanx,meany,meanz])

        #ranpoints[ranpoints>bs]-=bs
        #ranpoints[ranpoints<0.]+=bs


    if volumize==True:
        #########
        #VOLUMIZE
        #########
        f_array = np.zeros(np.shape(ranpoints))

        for i in range(nran):
            for face in voro[i]['faces']: 
                j = face['adjacent_cell']
                #print j
                if j<0: continue #negative adjacent cell means it's a wall
                dist = ranpoints[j]-ranpoints[i]
                
                dr = np.linalg.norm(dist)-msep

                r_ij = dist/np.linalg.norm(dist)
                f_array[i] += smooth*dr*r_ij 
            
            #voro = pyvoro.compute_voronoi(ranpoints,limits,msep)#,periodic=[True,True,True])
            
        ranpoints += f_array

    #ranpoints = copybox(ranpoints,buf,bs)
    voro = pyvoro.compute_voronoi(ranpoints,limits_cp,msep)
    
    
#TESTING: One last recenter
meanx=[]
meany=[]
meanz=[]
    
for i in range(len(ranpoints)):
   vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]

   meanx.append( np.mean(np.array(vert)[:,0]) )
   meany.append( np.mean(np.array(vert)[:,1]) )
   meanz.append( np.mean(np.array(vert)[:,2]) )

ranpoints = np.column_stack([meanx,meany,meanz])

#ranpoints[ranpoints>bs]-=bs
#ranpoints[ranpoints<0.]+=bs

#ranpoints = copybox(ranpoints,buf,bs)
voro = pyvoro.compute_voronoi(ranpoints,limits_cp,msep)


#VOLVEMOS AL BOX DE 1500
bs = 1500.
for axis in [0,1,2]:
    ranpoints = ranpoints[np.where(ranpoints[:,axis]<bs)]
    ranpoints = ranpoints[np.where(ranpoints[:,axis]>0.)]
limits = [[0.,bs],[0.,bs],[0.,bs]]

msep_real = (bs**nd/len(ranpoints))**(1./nd)
msep = 1.2*msep_real

##USO EL RBAR EN VEZ DEL MSEP!!!
rbar = (3.*bs**nd/(4.*np.pi*len(ranpoints)))**(1./3)
msep = 1.2*rbar

if periodicVoro==True: voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
if periodicVoro==False: voro = pyvoro.compute_voronoi(ranpoints,limits,msep)
    

vol2=[v.get('volume') for v in voro]

################
###
#Calculate Voronoi Volumes in CCVT for Comparison (MPE Cluster)
###

rcat = CSVCatalog('../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
rcat['Position']*=bs
real_mesh = rcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[bs,bs,bs],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
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
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkglass = r.power

gpoints = gcat['Position'].compute()
gvoro = pyvoro.compute_voronoi(gpoints,limits,msep,periodic=[True,True,True])
gvol=[v.get('volume') for v in gvoro]

################

plt.hist(vol1,bins=20,alpha=.7,label='Original',normed=False)
plt.hist(vol2,bins=20,alpha=.7,label='Recenter+Volumized',normed=False)
plt.hist(gvol,bins=20,alpha=.7,label='Glass',normed=False)
plt.hist(ccvt_vol,bins=20,alpha=.7,label='CCVT',normed=False)
plt.xlabel('Voronoi Cell Volume')
plt.legend()
plt.show()

points2 = ranpoints[np.where(ranpoints[:,2]<bs/10)]

plt.scatter(points2[:,0],points2[:,1],c='r',label='Recenter+Volumized')
plt.scatter(points1[:,0],points1[:,1],alpha=.4,label='Original')
plt.legend()
plt.show()

plt.scatter(ranpoints[:,0],ranpoints[:,1])
plt.show()

print 'Calculating P(k) - After'   
dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pk2 = r.power

#####
#Pks
#####
min = np.min(Pk2['power'].real)
print min
plt.hlines(min,Pk2['k'][0],Pk2['k'][-1],linestyle='-',label='Min = {}'.format(min))

plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
plt.loglog(Pk2['k'], Pk2['power'].real,label='After')# - Pk.attrs['shotnoise'])
plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass')
plt.loglog(Pk2['k'][:15],((bs**3/nran)*(Pk2['k'][:15]*msep/(2.*np.pi))**4),color='k',ls=':') 
plt.hlines((bs**nd)/nran, Pk1['k'][0],Pk2['k'][-1], colors='k',linestyles=':')
plt.vlines(2.*np.pi/msep,10,10E7,colors='k',linestyles=':',label='Mean Sep Voro')
plt.vlines(2.*np.pi/msep_real,10,10E7,colors='k',linestyles='-.',label='Mean Sep Real')
plt.legend()
plt.show()

#####################
#Ratio of Pk/Pk_ran
#####################
min = np.min(Pk2['power'].real/Pk1['power'].real)
print min
plt.hlines(min,Pk1['k'][0],Pk2['k'][-1],linestyle='-',label='Min = {}'.format(min))

#plt.loglog(Pk1['k'], Pk1['power'].real,label='Before')# - Pk.attrs['shotnoise'])
plt.loglog(Pk2['k'], Pk2['power'].real/Pk1['power'].real,label='After')# - Pk.attrs['shotnoise'])
plt.loglog(Pkccvt['k'], Pkccvt['power'].real/Pk1['power'].real,label='CCVT')
plt.loglog(Pkglass['k'], Pkglass['power'].real/Pk1['power'].real,label='Glass')
plt.hlines(1., Pk1['k'][0],Pk2['k'][-1], colors='k',linestyles=':')
plt.vlines(2.*np.pi/msep_real,10E-6,10,colors='k',linestyles=':',label='Mean Sep Real')
plt.vlines(2.*np.pi/msep,10E-6,10,colors='k',linestyles=':',label='Mean Sep Voro')
plt.vlines(2.*np.pi/bs,10E-6,10,colors='k',linestyles='-.',label='Box Size')

#plt.loglog(Pk2['k'],((bs**3/nran)*(Pk2['k']*msep/(2.*np.pi))**4),color='k',ls=':') 
plt.legend()
plt.show()


#Write Results
#results = {'min':[min], \
#           'recenter':[recenter], \
#           'volumize':[volumize], \
#           'smooth':[smooth], \
#           'niter':[niter], \
#           'pk':[Pk2['power'].real], \
#           'k':[Pk2['k']]} 
#results_t=Table(results)
#ascii.write(Table(results),'Results_{}.txt'.format(time.ctime()))
#f = open('Results_{}.txt'.format(time.ctime()),"w")
#f.write( str(results) )
#f.close()


#if volumize==True:
#    plt.plot(range(niter),norm,label=r'$\Delta r$')
#    plt.yscale('log')
#    plt.xscale('log')
#    plt.legend()
#    plt.show()

