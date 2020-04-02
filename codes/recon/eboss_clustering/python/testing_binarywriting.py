#tengo que correr el reconnew.py antes, de ahi saco el "new" que tiene las posiciones del catalogo zelrec

newbinary=bytearray(new)

newFile = open("prueba.bin", "wb")
newFile.write(newbinary)
newFile.close() 

dtype=[('Position', ('f4', 3))]
rcat1 = BinaryCatalog("prueba.bin",dtype)

rcat=rcat1.compute(rcat1['Position'])

plt.scatter(new[:,0],new[:,1],s=1)
plt.scatter(rcat[:,0],rcat[:,1],s=1)

##############
import struct

with open("prueba.bin", "wb") as file:
	file.write(struct.pack('{}d'.format(len(new)), *new))

##############
###
#ESTE FUNCIONA!
###
np.save("prueba.bin",new,allow_pickle=True)
                                                                   
table_ran=Table(np.load("prueba.bin.npy",allow_pickle=True),names=['x','y','z'])

rcat1 = ArrayCatalog(table_ran, BoxSize=[1,1,1])
rcat1['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])

rcat=rcat1.compute(rcat1['Position'])

plt.scatter(new[:,0],new[:,1],s=1)
plt.scatter(rcat[:,0],rcat[:,1],s=1)

