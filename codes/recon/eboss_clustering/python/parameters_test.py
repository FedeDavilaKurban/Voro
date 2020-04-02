plt.scatter(dat.x,dat.y)
plt.show()

plt.scatter(dat.newx,dat.newy)
plt.show()


for smooth_mult in [.5,1,2,3,5,8,10,20,40]:
    print(smooth_mult)
    exec(open("./reconnew.py").read())
    
    
for box in [1.,50.,100.,150.,200.,300.]:
    print(box)
    exec(open("./reconnew.py").read())
