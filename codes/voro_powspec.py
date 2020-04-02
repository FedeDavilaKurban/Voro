###
### 2D
###

nd = 2
nran = 100
bs = 1.
nbar = nran/bs**nd
msep = (bs**nd/nran)**(1./nd)

ranpoints = bs*np.random.uniform(size=(nran,nd))
limits = [[0.,bs],[0.,bs]]#,[0.,bs]]

voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])

#plt.scatter(ranpoints[:,0],ranpoints[:,1],c='g')
#for i in range(nran):
#    for x,y in voro[i]['vertices']:
#        plt.scatter(x,y)
#plt.show()

niter=100
for j in range(niter):
    print j
    if j==0: 
        plt.scatter(ranpoints[:,0],ranpoints[:,1])
        plt.show()
        
    voro = pyvoro.compute_2d_voronoi(ranpoints,limits,msep,periodic=[True,True])

    meanx=[]
    meany=[]
    for i in range(nran):
        vert=[(a[0],a[1]) for a in voro[i]['vertices']]
        meanx.append( np.mean(np.array(vert)[:,0]) )
        meany.append( np.mean(np.array(vert)[:,1]) )

       
    ranpoints = np.column_stack([meanx,meany])
    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    if j==(niter/2):
#        print ranpoints[ranpoints>bs]
#        print ranpoints[ranpoints<0.]
        plt.scatter(ranpoints[:,0],ranpoints[:,1])
        plt.show()

    
    if j==(niter-1):
#        print ranpoints[ranpoints>bs]
#        print ranpoints[ranpoints<0.]
        plt.scatter(ranpoints[:,0],ranpoints[:,1])
        plt.show()

#######################################################
###
### 3D
###

nd = 3
nran = 1000
bs = 1.
nbar = nran/bs**nd
msep = (bs**nd/nran)**(1./nd)

ranpoints = bs*np.random.uniform(size=(nran,nd))
limits = [[0.,bs],[0.,bs],[0.,bs]]

points=ranpoints[np.where(ranpoints[:,2]<.1)]
plt.scatter(points[:,0],points[:,1])
plt.show()

voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])


dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - Before'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.05)
Pk = r.power
plt.loglog(Pk['k'], Pk['power'].real,label='Random')# - Pk.attrs['shotnoise'])



niter=100
#mshift=[]
for j in range(niter):
    print j
#    if j==0:
#        points=ranpoints[np.where(ranpoints[:,2]<.2)]
#        plt.scatter(points[:,0],points[:,1])
#        plt.show()
        
    voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])

    meanx=[]
    meany=[]
    meanz=[]
    #shift=[]
    for i in range(nran):
        vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]
        meanx.append( np.mean(np.array(vert)[:,0]) )
        meany.append( np.mean(np.array(vert)[:,1]) )
        meanz.append( np.mean(np.array(vert)[:,2]) )
        #x,y,z = voro[i]['original']
        #shift.append( np.sqrt((x-meanx[-1])**2+(y-meany[-1])**2+(z-meanz[-1])**2) ) 

    #mshift.append( np.mean(np.array(shift)) )
    
    ranpoints = np.column_stack([meanx,meany,meanz])
    ranpoints[ranpoints>bs]-=bs
    ranpoints[ranpoints<0.]+=bs

    

#    if j==(niter/2):
#        points=ranpoints[np.where(ranpoints[:,2]<.2)]
#        plt.scatter(points[:,0],points[:,1])
#        plt.show()

    
#    if j==(niter-1):
#        points=ranpoints[np.where(ranpoints[:,2]<.2)]
#        plt.scatter(points[:,0],points[:,1])
#        plt.show()

#plt.plot(range(len(mshift)),mshift)
#plt.xlabel('N iteration')
#plt.ylabel('Mean Shift')
#plt.show()


dcat = ArrayCatalog(Table(ranpoints,names=['x','y','z']),BoxSize=bs)
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

print 'Calculating P(k) - After'    
real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=0.05)
Pk = r.power
plt.loglog(Pk['k'], Pk['power'].real,label='Voronized')# - Pk.attrs['shotnoise'])


plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
#plt.xlim(1E2, 1E3)
plt.legend()
plt.show()


points=ranpoints[np.where(ranpoints[:,2]<.1)]
plt.scatter(points[:,0],points[:,1])
plt.show()

