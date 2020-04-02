import numpy as np
from scipy.fftpack import fftfreq
import os
import pyfftw
import fastmodules
from nbodykit.lab import *
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.ndimage import gaussian_filter

class MiniCat:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.we = np.ones(len(x))
        self.size = x.size
        #self.dist = dist
        #self.x = x
        #self.y = y
        #self.z = z
        self.newx = x*1
        self.newy = y*1
        self.newz = z*1 

def get_shift(x, y, z, f_x, f_y, f_z, nbins, binsize):

    #xmin = 0.
    #ymin = 0.
    #zmin = 0.
    binsize1 = binsize
    nbins1 = nbins 
    
    xpos = (x-xmin)/binsize1
    ypos = (y-ymin)/binsize1
    zpos = (z-zmin)/binsize1
    
    i = xpos.astype(int)
    j = ypos.astype(int)
    k = zpos.astype(int)
    #print (k[np.where(k>=63)])
    ddx = xpos-i
    ddy = ypos-j
    ddz = zpos-k

    shift_x = np.zeros(x.size)
    shift_y = np.zeros(x.size)
    shift_z = np.zeros(x.size)

    for ii in range(2):
        for jj in range(2):
            for kk in range(2):

                weight = ( ((1-ddx)+ii*(-1+2*ddx))*\
                            ((1-ddy)+jj*(-1+2*ddy))*\
                            ((1-ddz)+kk*(-1+2*ddz))  )
                iii=i+ii
                jjj=j+jj
                kkk=k+kk
                #print(iii,jjj,kkk)
                
                #Para periodicidad
                iii[iii>=nbins1]-=nbins1
                jjj[jjj>=nbins1]-=nbins1
                kkk[kkk>=nbins1]-=nbins1
                
                pos = (iii, jjj, kkk)
                shift_x += f_x[pos]*weight
                shift_y += f_y[pos]*weight
                shift_z += f_z[pos]*weight

    return shift_x, shift_y, shift_z

niter=5
box = 100.
#box /= 2
inran = 16
nran = inran**3
np.random.seed(4251241242) 
data = np.random.uniform(low=-box,high=box,size=(nran,3))
rhom = nran/((box)**3)
msep = (((box)**3)/nran)**(1/3)

nbins = 128
smooth_mult = 1

binsize = box/nbins
smooth = binsize*smooth_mult
beta = 0.
bias = 1.
f = 0.
#box = bs
nthreads = 1
verbose = 0

write = 0
plot = 1 
compare = 0 #Comparar con los ccvt y glass de 16**3

dat = MiniCat(data[:,0], data[:,1], data[:,2])

dx = np.max(dat.x)-np.min(dat.x)
dy = np.max(dat.y)-np.min(dat.y)
dz = np.max(dat.z)-np.min(dat.z)
x0 = 0.5*(np.max(dat.x)+np.min(dat.x)) 
y0 = 0.5*(np.max(dat.y)+np.min(dat.y)) 
z0 = 0.5*(np.max(dat.z)+np.min(dat.z)) 

box = max([dx, dy, dz])+2*box
xmin = x0-box/2 
ymin = y0-box/2 
zmin = z0-box/2 
box = box
binsize = box/nbins
print('Box size [Mpc]:', box)
print('Nbins:', nbins)
print('Bin size [Mpc]:', binsize)
print('Bin size [1/Lbox]:', binsize/box)
print('Smooth_mult:', smooth_mult)
print('Num. of Particules:', nran)
print('Mean Part. Sep. [Mpc]:', msep)

plt.scatter(dat.x,dat.y)
plt.show()


for iloop in range(niter):
    print("Loop %d" % iloop)
    #-- Creating arrays for FFTW
    if iloop==0:
        delta  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        deltak = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        rho    = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        rhok   = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_x  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_y  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_z  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')

        #-- Initialize FFT objects and load wisdom if available
        wisdom_file = "wisdom."+str(nbins)+"."+str(nthreads)+'.npy'
        if os.path.isfile(wisdom_file) :
            print('Reading wisdom from ', wisdom_file)
            wisd = tuple(np.load(wisdom_file))
            print('Status of importing wisdom', pyfftw.import_wisdom(wisd))
        print('Creating FFTW objects...')
        fft_obj = pyfftw.FFTW(delta, delta, axes=[0, 1, 2], threads=nthreads)
        ifft_obj = pyfftw.FFTW(deltak, psi_x, axes=[0, 1, 2], \
                                threads=nthreads, \
                                direction='FFTW_BACKWARD')
        kr = fftfreq(nbins, d=binsize) * 2 * np.pi * smooth
        norm = np.exp(-0.5 * (  kr[:, None, None] ** 2 \
                                + kr[None, :, None] ** 2 \
                                + kr[None, None, :] ** 2))
    #else:
        #delta  = self.delta
        #deltak = self.deltak
        ##deltar = self.deltar
        #rho    = self.rho
        #rhok   = self.rhok
        #psi_x  = self.psi_x
        #psi_y  = self.psi_y
        #psi_z  = self.psi_z
        #fft_obj  = self.fft_obj
        #ifft_obj = self.ifft_obj
        #norm = self.norm
        
        
        
    #fft_obj = pyfftw.FFTW(delta, delta, threads=self.nthreads, axes=[0, 1, 2])
    #-- Allocate galaxies and randoms to grid with CIC method
    #-- using new positions
    if verbose:
        print('Allocating galaxies in cells...')
    deltag = np.zeros((nbins, nbins, nbins), dtype='float64')
    fastmodules.allocate_gal_cic(deltag, dat.newx, dat.newy, dat.newz, dat.we, dat.size, xmin, ymin, zmin, box, nbins, 1)
    if verbose:
        print('Smoothing...')
    #deltag = gaussian_filter(deltag, smooth/binsize)
    ##-- Smoothing via FFTs
    rho = deltag + 0.0j
    fft_obj(input_array=rho, output_array=rhok)
    fastmodules.mult_norm(rhok, rhok, norm)
    ifft_obj(input_array=rhok, output_array=rho)
    deltag = rho.real
    #print(deltag)
    if verbose:
        print('Computing density fluctuations, delta...')
    # normalize using the randoms, avoiding possible divide-by-zero errors
    #fastmodules.normalize_delta_survey(delta, deltag, self.alpha, self.ran_min)
    #del(deltag)  # deltag no longer required anywhere
    #rhom = len(dat.x)/box**3
    delta[:]  = deltag/rhom - 1.
    #w=np.where(deltar>self.ran_min)
    #delta[w] = delta[w]/(self.alpha*deltar[w])
    #w2=np.where((deltar<=self.ran_min)) 
    #delta[w2] = 0.
    #w3=np.where(delta>np.percentile(delta[w].ravel(), 99))
    #delta[w3] = 0.
    #del(w)
    #del(w2)
    #del(w3)
    del(deltag)
    if verbose:
        print('Fourier transforming delta field...')
    fft_obj(input_array=delta, output_array=delta)
    ## -- delta/k**2
    k = fftfreq(nbins, d=binsize) * 2 * np.pi
    fastmodules.divide_k2(delta, delta, k)
    if verbose:
        print('Inverse Fourier transforming to get psi...')
    fastmodules.mult_kx(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_x)
    fastmodules.mult_ky(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_y)
    fastmodules.mult_kz(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_z)

    # from grid values of Psi_est = IFFT[-i k delta(k)/(b k^2)], compute the values at the galaxy positions
    if verbose:
        print('Calculating shifts...')
    #print('Shape of psi_x:',np.shape(psi_x))
    shift_x, shift_y, shift_z = get_shift(dat.newx,dat.newy,dat.newz,psi_x.real,psi_y.real,psi_z.real,nbins,binsize)

    dat.newx += shift_x
    dat.newy += shift_y
    dat.newz += shift_z
    
    plt.scatter(dat.newx,dat.newy)
    plt.show()
    
    print(dat.newx[dat.newx>box])
    for axis in [dat.newx,dat.newy,dat.newz]:
        axis[axis>box/2] -= box/2
        axis[axis<-box/2] += box/2

    plt.scatter(dat.newx,dat.newy)
    plt.show()


    #if iloop%5==0: 
        #dcat = ArrayCatalog(Table(np.column_stack((dat.newx,dat.newy,dat.newz)),names=['x','y','z']),BoxSize=box)
        #dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
        #real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
        #r = FFTPower(real_mesh, mode='1d',dk=.07)
        #Pk = r.power
        #plt.loglog(Pk['k'], Pk['power'].real,label='{}'.format(iloop))


#for axis in [dat.x,dat.y,dat.z]:
#    axis *= box/4.
#    axis /= box/2.

############################################################################################
#plt.scatter(dat.newx,dat.newy)
#plt.show()

old = np.column_stack((dat.x,dat.y,dat.z))
new = np.column_stack((dat.newx,dat.newy,dat.newz))

dk = 7./box

dcat_old = ArrayCatalog(Table(old,names=['x','y','z']),BoxSize=box)
dcat_old['Position'] = transform.StackColumns(dcat_old['x'], dcat_old['y'], dcat_old['z'])
real_mesh_old = dcat_old.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r_old = FFTPower(real_mesh_old, mode='1d',dk=dk)
Pk_old = r_old.power

dcat_new = ArrayCatalog(Table(new,names=['x','y','z']),BoxSize=box)
dcat_new['Position'] = transform.StackColumns(dcat_new['x'], dcat_new['y'], dcat_new['z'])
real_mesh_new = dcat_new.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
r_new = FFTPower(real_mesh_new, mode='1d',dk=dk)
Pk_new = r_new.power

if compare:
    if Pkccvt is None:
        rcat = CSVCatalog('../../../../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
        rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
        rcat['Position']*=box
        real_mesh = rcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
        r = FFTPower(real_mesh, mode='1d',dk=dk)
        Pkccvt = r.power

    if Pkglass is None:
        gcat = CSVCatalog('../../../../data/glass_16.dat',names=['x','y','z'])
        gcat['Position'] = transform.StackColumns(gcat['x'], gcat['y'], gcat['z'])
        gcat['Position']*=box
        real_mesh = gcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
        r = FFTPower(real_mesh, mode='1d',dk=dk)
        Pkglass = r.power


if plot:
    plt.loglog(Pk_old['k'], Pk_old['power'].real,label='Random')
    plt.loglog(Pk_new['k'], Pk_new['power'].real,label='Zel. Recon.')
    if compare:
        plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
        plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass')
    #plt.loglog(Pk_old['k'][:15],((box**3/nran)*(Pk_old['k'][:15]*msep/(2.*np.pi))**4),color='k',ls=':') 
    plt.legend()
    plt.show()

############################################################################################
if write:
    print('Writing')
    ascii.write(Table(old,names=['x','y','z']),output='../../../../data/newCat_{}_{}_{}_{}_{}.txt'.format(niter,nbins,box,inran,smooth_mult),overwrite=True)
    ascii.write(Table(np.column_stack((Pk_new['k'].real,Pk_new['power'].real)),names=['k','Pk']),output='../../../../data/newPk_{}_{}_{}_{}_{}.txt'.format(niter,nbins,box,inran,smooth_mult),overwrite=True)
