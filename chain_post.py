import matplotlib.pyplot as plt
import os
import numpy as np

class PostPlt:
    def __init__(self,chainrst='test',xdatanm=['section','xx'],ydatanm=['section','yy']):
        self.chainrst = chainrst
        self.xdatanm = xdatanm
        self.ydatanm = ydatanm
        self.blocksnm = []
        self.xdataar = []
        self.ydataar = []
        for dir in os.listdir(chainrst):
            if dir[:6] == 'chain_':
                self.blocksnm.append(dir)
        for blocknm in self.blocksnm:

            self.xdataar.append(np.loadtxt(chainrst+'/'+blocknm+'/'+xdatanm[0]+'/'+xdatanm[1]))
            self.ydataar.append(np.loadtxt(chainrst+'/'+blocknm+'/'+ydatanm[0]+'/'+ydatanm[1]))
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(1,1,1)

    def adline(self,xdata,ydata,ax):
        ln, = ax.plot(xdata,ydata)
        return ln

    def plotprocess(self):
        self.lines = []
        for xdata, ydata in zip(self.xdataar,self.ydataar):
            self.lines.append(self.adline(xdata,ydata,self.axes))
            self.lines[-1].set_label(self.blocksnm[len(self.lines)-1].replace(self.xdatanm[0],'').replace(self.ydatanm[0],''))
        self.axes.legend(loc=2)
        self.fig.savefig(self.chainrst+'/figure.png')


def main():
    pltins = PostPlt('testdemo1',['distances','z.txt'],['distances','a.txt'])
    pltins.plotprocess()
    print pltins.xdataar, pltins.ydataar

if __name__ == "__main__":
    main()
