from pylab import plot,show,fill,legend,hlines,title,xlabel,ylabel
def par(p,s):
    return 1.0/((1.0-p)+p/s)
dashes = ['--','-.',':','-']
plot(range(1,1000),  [par(.57835,s) for s in range(1,1000)],label="Our Theoretical Speedup (57.835%)")
hlines(par(.57835,41),0,1000,label="Theoretical Speedup")
k = 0
for i in range(50,100,20):
    plot(range(1,1000),[par(i/100.0,s) for s in range(1,1000)],dashes[k],label=str(i)+ "%")
    k=k+1
plot(range(1,1000),[par(.95,s) for s in range(1,1000)],'--',label=str(95)+ "%")
#plot(range(1,1000),[par(.97,s) for s in range(1,1000)],'--',label=str(100-97)+ "% Serial Code")
fill([40,40,65,65], [0,40,40,0], 'b', alpha=0.2, edgecolor='r')

from matplotlib.font_manager import fontManager, FontProperties
font= FontProperties(size='x-small');
legend(title="Parallelism",loc="lower right")
title("Total Speedup as a function of Parallel Speedup")
ylabel("Total Speedup")
xlabel("Parallel Speedup")
show()
