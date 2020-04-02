"""
Test

Averiguando que pasa si genero un box mas grande de lo que necesito para evitar problemas de borde

Este codigo no tiene sentido ya, despues de haber creado la funcion "copybox" con el codigo del Enrique. Es mas correcto de esa forma

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import pyvoro
from nbodykit.lab import *
import time

recenter = True
volumize = True
periodicVoro = False
smooth = .2 #se usa en el volumizado 
niter=10

nd = 3
nran = 16**3
bs = 1500.

#######
#TEST generando un box mayor para evitar problemas de borde
rho = nran/bs**3
bs = 1700.
nran = int(rho*bs**3)
#######


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
#ranpoints = -100.+(bs)*np.random.uniform(size=(nran,nd)) #TESTING


limits = [[0.,bs],[0.,bs],[0.,bs]]
#limits = [[-100.,bs],[-100.,bs],[-100.,bs]] #TESTING

points1 = ranpoints[np.where(ranpoints[:,2]<bs/10)]

if periodicVoro==True: voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
if periodicVoro==False: voro = pyvoro.compute_voronoi(ranpoints,limits,msep)

dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - Before'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pk1 = r.power

#t1 = time.time()

vol1 = [v.get('volume') for v in voro]
norm = np.zeros(niter)


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

        for i in range(nran):
            vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]

            meanx.append( np.mean(np.array(vert)[:,0]) )
            meany.append( np.mean(np.array(vert)[:,1]) )
            meanz.append( np.mean(np.array(vert)[:,2]) )

        ranpoints = np.column_stack([meanx,meany,meanz])

        ranpoints[ranpoints>bs]-=bs
        ranpoints[ranpoints<0.]+=bs

    if periodicVoro==True: voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
    if periodicVoro==False: voro = pyvoro.compute_voronoi(ranpoints,limits,msep)

    if volumize==True:
        #########
        #VOLUMIZE
        #########
        f_array = np.zeros(np.shape(ranpoints))
        #f_array = np.zeros((np.shape(ranpoints[0])))
        #print 'Initial f_array:',f_array
        #dr_array = np.zeros(np.shape(ranpoints))
        for i in range(nran):
            for face in voro[i]['faces']: 
                j = face['adjacent_cell']
                #print j
                if j<0: continue #negative adjacent cell means it's a wall
                #if j==i: continue
                dist = ranpoints[j]-ranpoints[i]
                
                dr = np.linalg.norm(dist)-msep
                                
                #if dr<0.: neg+=1
                #if dr>0.: pos+=1
                #print dr
                r_ij = dist/np.linalg.norm(dist)
                f_array[i] += smooth*dr*r_ij #lo reemplace por una sola f, para cambiar la posicion una por una
                #print 'f:',smooth*dr*r_ij
                #f_array += smooth*dr*r_ij
                #print 'f_array:',f_array
                #f = smooth*dr*r_ij
                #print i,idx,dist,r_ij,dr_array[idx]
            
            #print f_array[i]
            #ranpoints[i] += f_array
            #ranpoints[i] += f
            #if np.any(np.isnan(ranpoints[i]))==True:
            #    print 'nan',i
            
            #voro = pyvoro.compute_voronoi(ranpoints,limits,msep)#,periodic=[True,True,True])
            
        ranpoints += f_array
        #dr_array = f_array
        #norm[_] = np.mean(np.linalg.norm(dr_array,axis=1))

        #ranpoints += dr_array #lo reemplaze por el 'ranpoints[i] += f' mas arriba

        #Check if there are particles outside of box
        #print 'Percentage of Out Of Box Particles:',(len(ranpoints[ranpoints>bs])+len(ranpoints[ranpoints<0.]))*100/len(ranpoints),'%'
        ranpoints[ranpoints>bs]-=bs
        ranpoints[ranpoints<0.]+=bs
        #print 'Percentage of Out Of Box Particles:',(len(ranpoints[ranpoints>bs])+len(ranpoints[ranpoints<0.]))*100/len(ranpoints),'%'

    
    if periodicVoro==True: voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
    if periodicVoro==False: voro = pyvoro.compute_voronoi(ranpoints,limits,msep)

#TESTING: One last recenter #Works. Not testing
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

#TESTING Volvemos al box de 1500###############
bs = 1500
for axis in [0,1,2]:
    ranpoints = ranpoints[np.where(ranpoints[:,axis]<bs)]
    ranpoints = ranpoints[np.where(ranpoints[:,axis]>100.)]
ranpoints -= 100
limits = [[0.,bs],[0.,bs],[0.,bs]]

msep_real = (bs**nd/len(ranpoints))**(1./nd)
msep = 1.2*msep_real
if periodicVoro==True: voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])
if periodicVoro==False: voro = pyvoro.compute_voronoi(ranpoints,limits,msep)
#############################

#t2 = time.time()

#ascii.write(ranpoints,'../volumized_cat/volumized_{}_{}.txt'.format(nran,niter),format='no_header',overwrite=True)

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
