nran=100
np.random.seed(2022342)
ranpoints = bs*np.random.uniform(size=(nran,nd))


plt.scatter(ranpoints[:,0],ranpoints[:,1])
plt.show()

pos = ranpoints
buf=150.
box=bs
#pos = np.genfromtxt(fin)[:, :3]

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

plt.scatter(pos[:,0],pos[:,1])
plt.show()
