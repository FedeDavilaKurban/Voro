import os
import sys
import numpy as np
from nbodykit.lab import *

from astropy.io import ascii
from astropy import table
from astropy.table import Table

def get_xi0246(corr,nbins_m,nbins_s):

    xi_sm = corr.corr.data['corr']

    ##Modificado para que lea los valores de entrada
    ##nbins_m=30 # number of bins in mu
    ##nbins_s=29 # number of bins in s
    dmu=1.0/nbins_m

    rs = corr.D1D2.coords['r']
    mu = corr.D1D2.coords['mu']

    xi_s0 = np.zeros(nbins_s)
    xi_s2 = np.zeros(nbins_s)
    xi_s4 = np.zeros(nbins_s)
    xi_s6 = np.zeros(nbins_s)

    sr = np.zeros(nbins_s)
    rm = np.zeros(nbins_m)

    l0 = 0.0
    l1 = 1.0
    l2 = 2.0
    l3 = 3.0

    for i in range(nbins_s):
        sr[i] = rs[i]
        for j in range(nbins_m):
                rm[j]=mu[j]
                xi_s0[i]  += (4.0*l0+1.0)*xi_sm[i,j]*1.0*dmu
                xi_s2[i]  += (4.0*l1+1.0)*xi_sm[i,j]*((3*rm[j]**2 - 1.0)/2.0)*dmu
                xi_s4[i]  += (4.0*l2+1.0)*xi_sm[i,j]*((35*rm[j]**4 - 30*rm[j]**2 + 3.0)/8.0)*dmu
                xi_s6[i]  += (4.0*l3+1.0)*xi_sm[i,j]*((231*rm[j]**6 - 315*rm[j]**4 + 105*rm[j]**2 - 5)/16.0)*dmu

    return xi_s0, xi_s2, xi_s4, xi_s6
##################################################################################
i = int(os.environ['SGE_TASK_ID'])*2-2
print i

bs = 1500.
nbins_m = 30
nbins_s = 29
names=['xi_0','xi_2','xi_4','xi_6']

os.chdir('../data')

datafile = 'Galaxies_HOD_001_z0.57.dat'

dcat = CSVCatalog(datafile,names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

###Defino el dtype para leer el archivo binario
###Al definirlo asi ya le doy el atributo 'Position' para que despues calcule la correlacion
dtype=[('Position', ('f4', 3))]

if int(sys.argv[1]) == 1: zelrec_dir = 'zelrec/1x87/'
if int(sys.argv[1]) == 2: zelrec_dir = 'zelrec/2x87/'
if int(sys.argv[1]) == 4: zelrec_dir = 'zelrec/4x87/'
if int(sys.argv[1]) == 8: zelrec_dir = 'zelrec/8x87/'

#files = sorted(os.listdir(zelrec_dir))
files = sorted(os.listdir(zelrec_dir))
os.chdir(zelrec_dir)

f1 = files[i]
f2 = files[i+1]

print f1,f2

###Calculo y escribo los xi_l
#id1=f1[-7:-4]
#id2=f2[-7:-4]
#print '~~~~~~~~~~~~~~~~~~~~~~~~',id1,id2,'~~~~~~~~~~~~~~~~~~~~~~~~~'

zr1 = ascii.read(f1,names=['x','y','z'])
zr2 = ascii.read(f2,names=['x','y','z'])

for col in ['x','y','z']:
    zr1[col]*=bs
    zr2[col]*=bs

rcat1 = ArrayCatalog(zr1, BoxSize=[bs,bs,bs])
rcat1['Position'] = transform.StackColumns(rcat1['x'], rcat1['y'], rcat1['z'])
nran = rcat1.size

rcat2 = ArrayCatalog(zr2, BoxSize=[bs,bs,bs])
rcat2['Position'] = transform.StackColumns(rcat2['x'], rcat2['y'], rcat2['z'])


corr = SimulationBox2PCF('2d',data1=dcat,data2=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat1,randoms2=rcat2,BoxSize=[bs,bs,bs], periodic=True)

xi_l = get_xi0246(corr,nbins_m,nbins_s)

os.chdir('../../out/zelrec/')

wfile = 'xi_l_zelrec_{}_{}-{}.txt'.format(rcat1.size,i,i+1)

ascii.write(xi_l,wfile,overwrite=True,format='no_header')
