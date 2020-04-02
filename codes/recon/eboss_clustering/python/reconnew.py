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
import time
import seaborn

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

    xmin = 0.
    ymin = 0.
    zmin = 0.
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

niter=100
box = 1.
inran = 87
nran = inran**3
seed=np.random.seed(lll) 
data = box*np.random.uniform(size=(nran,3))

nbins = 696
smooth_mult = 1

binsize = box/nbins
smooth = binsize*smooth_mult

rhom_t = nran/((box)**3)
rhom = rhom_t*(binsize)**3
msep = rhom_t**(-1/3) #celda
msep = box/((nran*4*np.pi/3)**(1/3))

#testing
smooth = msep*.5
#smooth = binsize
#testing

print('smooth=',smooth,'smoothbefore=',binsize*smooth_mult)


beta = 0.
bias = 1.
f = 0.
nthreads = 1
verbose = 0

write = 0
plot = 1
compare = 1 #Comparar con los ccvt y glass de 16**3

dat = MiniCat(data[:,0], data[:,1], data[:,2])

xmin=np.min(dat.x)
ymin=np.min(dat.y)
zmin=np.min(dat.z)
#print('Box size [Mpc]:', box)
#print('Nbins:', nbins)
#print('Bin size [Mpc]:', binsize)
##print('Bin size [1/Lbox]:', binsize/box)
#print('Smooth_mult:', smooth_mult)
#print('Num. of Particules:', nran)
##print('Mean Part. Sep. [Mpc]:', msep)


#For some plots
c=seaborn.light_palette('orange',10)
iplot=0
dk = 7./box


#plt.scatter(dat.x[np.where(dat.newz<box/10)],dat.y[np.where(dat.newz<box/10)])#,color=c[-1])
#plt.xlim((-.05,1.05))
#plt.ylim((-.05,1.05))
##plt.savefig('../../../../plots/cat_forgif_{}.png'.format(iplot))
##plt.close()
#plt.show()


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

    if np.all(psi_x.real[0] == np.zeros(np.shape(psi_x.real[0]))): 
        print('No displacement field. Breaking.')
        break


    # from grid values of Psi_est = IFFT[-i k delta(k)/(b k^2)], compute the values at the galaxy positions
    if verbose:
        print('Calculating shifts...')
    #print('Shape of psi_x:',np.shape(psi_x))
    shift_x, shift_y, shift_z = get_shift(dat.newx,dat.newy,dat.newz,psi_x.real,psi_y.real,psi_z.real,nbins,binsize)
    
    dat.newx += shift_x
    dat.newy += shift_y
    dat.newz += shift_z
    
    #iplot+=1
    #plt.scatter(dat.newx[np.where(dat.newz<box/10)],dat.newy[np.where(dat.newz<box/10)])#,color=c[iloop+1])
    #plt.xlim((-.05,1.05))
    #plt.ylim((-.05,1.05))
    #plt.savefig('../../../../plots/cat_forgif_{}.png'.format(iplot))
    #plt.close()
    ##plt.show()

    #print(dat.newx[dat.newx>2*box])
    for axis in [dat.newx, dat.newy, dat.newz]:
        axis[axis>box] -= box
        axis[axis<0.] += box
    
    #iplot+=1
    #plt.scatter(dat.newx[np.where(dat.newz<box/10)],dat.newy[np.where(dat.newz<box/10)])#,color=c[iloop+1])
    #plt.xlim((-.05,1.05))
    #plt.ylim((-.05,1.05))
    #plt.savefig('../../../../plots/cat_forgif_{}.png'.format(iplot))
    #plt.close()
    ##plt.show()

    
    #if iloop%1==0: 
        #dcat = ArrayCatalog(Table(np.column_stack((dat.newx,dat.newy,dat.newz)),names=['x','y','z']),BoxSize=box)
        #dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
        #real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256)
        #r = FFTPower(real_mesh, mode='1d',dk=dk)
        #Pk = r.power
        #plt.loglog(Pk['k'], Pk['power'].real,color=c[iplot],label='{}'.format(iloop))
        #iplot += 1

    wisdom_file = "wisdom."+str(nbins)+"."+str(nthreads)+'.npy'
    if iloop==0 and not os.path.isfile(wisdom_file):
        wisd=pyfftw.export_wisdom()
        np.save(wisdom_file, wisd)
        print('Wisdom saved at', wisdom_file)

############################################################################################
#plt.scatter(dat.newx,dat.newy)
#plt.show()

old = np.column_stack((dat.x,dat.y,dat.z))
new = np.column_stack((dat.newx,dat.newy,dat.newz))


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
	if inran==16:
		#if 'Pkccvt' not in locals():
			rcat = CSVCatalog('../../../../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
			rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
			rcat['Position']*=box
			real_mesh = rcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
			r = FFTPower(real_mesh, mode='1d',dk=dk)
			Pkccvt = r.power

		#if 'Pkglass' not in locals():
			gcat = CSVCatalog('../../../../data/glass_16.dat',names=['x','y','z'])
			gcat['Position'] = transform.StackColumns(gcat['x'], gcat['y'], gcat['z'])
			gcat['Position']*=box
			real_mesh = gcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
			r = FFTPower(real_mesh, mode='1d',dk=dk)
			Pkglass = r.power
		    
	if inran==87:
		#if 'Pkglass' not in locals():
			dtype=[('Position', ('f4', 3))]
			gcat=BinaryCatalog('../../../../data/glass_001.bin',dtype,header_size=4)
			gcat['Position']/=1500.
			real_mesh = gcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
			r = FFTPower(real_mesh, mode='1d',dk=dk)
			Pkglass = r.power
	
		#if 'Pkccvt' not in locals():
			rcat = CSVCatalog('../../../../data/ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
			rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
			rcat['Position']*=box
			real_mesh = rcat.to_mesh(compensated=True, window='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
			r = FFTPower(real_mesh, mode='1d',dk=dk)
			Pkccvt = r.power

if plot:
    plt.loglog(Pk_old['k'][:40],((box**3/nran)*(Pk_old['k'][:40]*msep*.4/(2.*np.pi))**4),color='k',ls=':') 
    plt.loglog(Pk_old['k'], Pk_old['power'].real,label='Random')
    if compare:
        plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
        plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass')
    plt.loglog(Pk_new['k'], Pk_new['power'].real,label='Zel. Recon.')
    #plt.loglog(Pk_old['k'][:15],((box**3/nran)*(Pk_old['k'][:15]*msep/(2.*np.pi))**4),color='k',ls=':')
    plt.legend()
    #plt.ylim((8*10E-13,2*10E-4))
    plt.savefig('../../../../plots/pk_{}_{}_{}_{}.png'.format(inran,niter,nbins,smooth_mult))
    #plt.close()
    plt.show()

############################################################################################
if write:
	filename = "zelrec_{0:0=3d}".format(lll) #Ese formato es para escribir 3 digitos
	print('Writing {}'.format(filename))
    #ascii.write(Table(old,names=['x','y','z']),output='../../../../data/ZelRec_{}_{}_{}_{}.txt'.format(niter,nbins,inran,smooth_mult),overwrite=True)
    #ascii.write(Table(np.column_stack((Pk_new['k'].real,Pk_new['power'].real)),names=['k','Pk']),output='../../../../data/ZelRecPk_{}_{}_{}_{}.txt'.format(niter,nbins,inran,smooth_mult),overwrite=True)
    ###
    #ESTE DE ABAJO ES EL QUE ESTOY USANDO
    ###
    #ascii.write(Table(old,names=['x','y','z']),output='../../../../data/zelrec/ZelRec_{}_{}_{}_{}_{}.txt'.format(niter,nbins,inran,smooth_mult,seed,lll),overwrite=True)
    #ascii.write(Table(np.column_stack((Pk_new['k'].real,Pk_new['power'].real)),names=['k','Pk']),output='../../../../data/zelrec/ZelRecPk_{}_{}_{}_{}_{}.txt'.format(niter,nbins,inran,smooth_mult,lll),overwrite=True)
    
    
	np.save("../../../../data/zelrec/{}".format(filename),new,allow_pickle=True) 

