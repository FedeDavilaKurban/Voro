files = !ls ../../../../data/*Pk*256*87*

for f in files:
    Pk = ascii.read(f,names=['k','Pk'],)
    plt.loglog(Pk['k'], Pk['Pk'])
plt.show()
