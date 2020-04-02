"""
Ver como varia el mean shift con la cantidad de particulas.

Cuantas iteraciones hace falta para que el mean shift caiga a 0
"""

nd = 3
bs = 1.
limits = [[0.,bs],[0.,bs],[0.,bs]]

niter=100


n_mshift=[]
n_mshift_msep=[]
for nran in [100,1000,5000]:
    print 'Nran=',nran
    nbar = nran/bs**nd
    msep = (bs**nd/nran)**(1./nd)

    ranpoints = bs*np.random.uniform(size=(nran,nd))

    voro = pyvoro.compute_voronoi(ranpoints,limits,msep,periodic=[True,True,True])


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
        shift=[]
        for i in range(nran):
            vert=[(a[0],a[1],a[2]) for a in voro[i]['vertices']]
            meanx.append( np.mean(np.array(vert)[:,0]) )
            meany.append( np.mean(np.array(vert)[:,1]) )
            meanz.append( np.mean(np.array(vert)[:,2]) )
            x,y,z = voro[i]['original']
            shift.append( np.sqrt((x-meanx[-1])**2+(y-meany[-1])**2+(z-meanz[-1])**2) ) 

        n_mshift.append( np.mean(np.array(shift)) )
        n_mshift_msep.append( np.mean(np.array(shift))/msep )
        
        ranpoints = np.column_stack([meanx,meany,meanz])
        ranpoints[ranpoints>bs]-=bs
        ranpoints[ranpoints<0.]+=bs


plt.plot(range(niter),n_mshift[:niter],label='Nran=100')
plt.plot(range(niter),n_mshift[niter:2*niter],label='Nran=1000')
plt.plot(range(niter),n_mshift[2*niter:3*niter],label='Nran=5000')
plt.legend()
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('N iteration')
plt.ylabel('Mean Shift')
plt.show()

plt.plot(range(niter),n_mshift_msep[:niter],label='Nran=100')
plt.plot(range(niter),n_mshift_msep[niter:2*niter],label='Nran=1000')
plt.plot(range(niter),n_mshift_msep[2*niter:3*niter],label='Nran=5000')
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.xlabel('N iteration')
plt.ylabel('Mean Shift / Mean Separation')
plt.show()
