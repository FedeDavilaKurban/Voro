"""
Test

Intentando optimizar el calculo del volumizado.
La idea es no calcular las distancias vecino a vecino, si no calcular de una las distancias a todos los vecinos.

Funciona

P.D.: ?Quizas se podria hacer un paso mas y hacerlo para todas las particulas de una, en vez de iterar particula por particula?
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import pyvoro
from nbodykit.lab import *
import time

smooth = .5 #se usa en el volumizado 
niter=10
buf=150.

nd = 3
nran = 16**3
bs = 1500.

np.random.seed(2022342)
ranpoints = bs*np.random.uniform(size=(nran,nd))

##USO EL RBAR EN VEZ DEL MSEP!!!
rbar = (3.*bs**nd/(4.*np.pi*len(ranpoints)))**(1./3)
msep = 1.2*rbar

limits = [[0.,bs],[0.,bs],[0.,bs]]

voro = pyvoro.compute_voronoi(ranpoints,limits,msep)


f_array = np.zeros(np.shape(ranpoints))

t1=time.time()
for i in range(len(ranpoints)):
    for face in voro[i]['faces']: 
        j = face['adjacent_cell']
        #print j
        if j<0: continue #negative adjacent cell means it's a wall
        dist = ranpoints[j]-ranpoints[i]
        
        dr = np.linalg.norm(dist)-msep

        r_ij = dist/np.linalg.norm(dist)
        f_array[i] += smooth*dr*r_ij 
t2=time.time()
print f_array[0]

f_array = np.zeros(np.shape(ranpoints))
f1 = f_array
t3=time.time()
for i in range(len(ranpoints)):
    
    js = np.array([ii['adjacent_cell'] for ii in voro[i]['faces']] )
    js = js[js>0]
    dists = ranpoints[js]-ranpoints[i]
    norms = np.linalg.norm(dists,axis=1)
    drs = norms-msep

    r_ijs = dists/norms[:,None]
    f_array[i] = np.sum(smooth*drs[:,None]*r_ijs,axis=0)
t4=time.time()
    
print f_array[0]
print t2-t1,t4-t3
