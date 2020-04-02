"""
Volumizando usando el Voronoi de Scipy

16/09/2019 lo toqueteé para el Seminario 1
"""
##########
def copybox2d(pos,buf,box):

    g1 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box])
    g1[:, 0] -= box
    g1[:, 1] -= box
    #g1[:, 2] -= box
    g2 = np.asarray([i for i in pos if box-buf<i[1]<box])
    g2[:, 1] -= box
    #g2[:, 2] -= box
    g3 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box])
    g3[:, 0] += box
    g3[:, 1] -= box
    #g3[:, 2] -= box
    g4 = np.asarray([i for i in pos if box-buf<i[0]<box])
    g4[:, 0] -= box
    #g4[:, 2] -= box
    g5 = np.asarray([i for i in pos if 0<i[1]<buf])
    g5[:, 1] += box
    g6 = np.asarray([i for i in pos if 0<i[0]<buf])
    g6[:, 0] += box
    #g6[:, 2] -= box
    g7 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf])
    g7[:, 0] -= box
    g7[:, 1] += box
    #g7[:, 2] -= box
    g8 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf])
    g8[:, 0] += box
    g8[:, 1] += box


    pos = np.vstack([pos, g1, g2, g3, g5, g4, g6, g7,g8])
    return pos
###########


from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import matplotlib.pyplot as plt

bs=10.
buf=3.
np.random.seed(2022342)
ranpoints=bs*np.random.uniform(size=(100,2))
plt.scatter(ranpoints[:,0],ranpoints[:,1])
plt.savefig("voro2d_0.png")
#plt.show()

ranpoints = copybox2d(ranpoints,buf,bs)
msep = (100./bs**3)**(-1./3)
smooth = 0.2
#plt.scatter(ranpoints[:,0],ranpoints[:,1])
#plt.show()

vor = Voronoi(ranpoints)

fig = voronoi_plot_2d(vor)
plt.xlim(0.,bs)
plt.ylim(0.,bs)
plt.savefig("voro2d_1.png")
#plt.show()

niter = 10

for _  in range(niter):
    print _
    f_array = np.zeros(np.shape(ranpoints))
    for i in range(len(ranpoints)):
        
        pairs=vor.ridge_points[np.where(vor.ridge_points==i)[0]]
        js = pairs[pairs!=i]
        dists = ranpoints[js]-ranpoints[i]
        norms = np.linalg.norm(dists,axis=1)
        drs = norms-msep

        r_ijs = dists/norms[:,None]
        f_array[i] = np.sum(smooth*drs[:,None]*r_ijs,axis=0)

    ranpoints += f_array
    ranpoints[ranpoints>(bs+buf)]-=(bs+buf)
    ranpoints[ranpoints<-buf]+=(bs+buf)
    #plt.scatter(ranpoints[:,0],ranpoints[:,1])
    #plt.show()

    vor = Voronoi(ranpoints)
    fig = voronoi_plot_2d(vor)
    plt.xlim(0.,bs)
    plt.ylim(0.,bs)
    plt.savefig("voro2d_{}.png".format(_+2))
    #plt.show()

_+=3
print 'Back to original box size'
for axis in [0,1]:
    ranpoints = ranpoints[np.where(ranpoints[:,axis]<bs)]
    ranpoints = ranpoints[np.where(ranpoints[:,axis]>0.)]
vor = Voronoi(ranpoints)
fig = voronoi_plot_2d(vor)
plt.xlim(0.,bs)
plt.ylim(0.,bs)
plt.savefig("voro2d_{}.png".format(_))
#plt.show()


#dcat = ArrayCatalog(Table(ranpoints,names=['x','y']),BoxSize=bs)
#dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'])
#real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
#r = ProjectedFFTPower(real_mesh, dk=dk, axes=[0,1])
#Pk2 = r.power
#plt.loglog(Pk2['k'], Pk2['power'].real)#,label='{}'.format(_))
#plt.show()
