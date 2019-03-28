#!/usr/bin/python
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import numpy
import numpy as np
import os
import fnmatch

VAR1=os.getenv('VAR1','Cost Function')
VAR2=os.getenv('VAR2','Gradient Norm Reduction')
file_name = 'costgrad.txt'  # output from grep and paste command

def integers(a, b):
    return list(range(a, b+1))

def readdata():

    dalogfiles = []
    for files in os.listdir('../testoutput/'):
        if fnmatch.fnmatch(files, '?d*.test.log.out'):
            dalogfiles.append('../testoutput/'+files)
            #print(dalogfiles)

    for dalogfile in dalogfiles:
        print(dalogfile)
        file_name = 'costgrad_'+dalogfile[14:][:-13]+'.txt'

        exists = os.path.isfile(file_name)
        if exists:
            os.system('rm '+ file_name)

        #grep iterations and cost function.
        cmd = 'grep "Quadratic cost function: J " '+ dalogfile +' \
            | grep -o -P "(\(\K[^\)]+)|(=\K.+)"  |paste -d " " - - > cost.txt'
        os.system(cmd)

        #grep "Norm reduction" 
        cmd = 'grep "DRIPCG end of iteration" '+ dalogfile +' \
            | grep -o -P "=\K.+" > grad.txt'
        os.system(cmd)

        #combine cost and grad:
        cmd = 'paste cost.txt grad.txt > '+ file_name +' ; rm cost.txt grad.txt' 
        os.system(cmd)

        alist = open(file_name).read().split()
        iters  = numpy.asarray(alist[0::3])  
        cost = numpy.asarray(alist[1::3])
        grad = numpy.asarray(alist[2::3])
        forx=integers(1,len(iters)) 
        plot(forx,cost,iters,VAR1,dalogfile[14:])
        plot(forx,grad,iters,VAR2,dalogfile[14:])

def plot(forx,value,iters,VAR,expname): 
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Iterations',fontsize=16)
    ax1.set_xticks(forx[::2])
    ax1.set_xticklabels(iters.astype(np.int)[::2]) #,rotation=45)
    ax1.set_ylabel('%s'%VAR,fontsize=16)

    plt.plot(forx,value,'r.-')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,4))
    plt.grid(True)
    plt.savefig('%s_%s.png'%(VAR[:4],expname[:-13]),dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
