#!/usr/bin/python
import matplotlib
matplotlib.use('AGG')
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import numpy
import numpy as np
import os
import fnmatch

'''
Directory Structure:
test/
 ├── testoutput/
 │   ├── 3dvar.run
 │   ├── 3denvar_2stream_bumploc_unsinterp.run
 │   ├── ...
 ├── graphics/
 │   ├── plot_cost_grad.py
 │   ├── ...
'''

VAR1=os.getenv('VAR1','Quadratic cost function')
#VAR2=os.getenv('VAR2','Gradient reduction')
VAR2=os.getenv('VAR2','Norm reduction')
file_name = 'costgrad.txt'  # output from grep and paste command

def integers(a, b):
    return list(range(a, b+1))

def readdata():

    dalogfiles = []
    for files in os.listdir('../testoutput/'):
        if fnmatch.fnmatch(files, '?d*.run'):
            dalogfiles.append('../testoutput/'+files)
            #print(dalogfiles)

    for dalogfile in dalogfiles:
        print('check dalogfile=',dalogfile,dalogfile[14:][:-4])
        file_name = 'costgrad_'+dalogfile[14:][:-4]+'.txt'

        exists = os.path.isfile(file_name)
        if exists:
            os.system('rm '+ file_name)

        #grep iterations and cost function.
        cmd = 'grep "'+VAR1+': J " '+ dalogfile +' \
            | grep -o -P "(\(\K[^\)]+)|(=\K.+)"  |paste -d " " - - > cost.txt'
        os.system(cmd)

        cmd = 'grep "'+VAR2+'" '+ dalogfile +' \
            | grep -o -P "=\K.+" > gradorig.txt'
        os.system(cmd)

        nLineCost = int(os.popen('wc -l < cost.txt').read())
        nLineGrad = int(os.popen('wc -l < gradorig.txt').read())

        #some outputs included GMRESR, remove 'Gradient reduction' or 'Norm reduction' from it.
        if nLineCost !=  nLineGrad:
            cmd = 'head -n '+str(nLineCost)+' gradorig.txt > grad.txt'
        else:
            cmd = 'cp gradorig.txt grad.txt'
        os.system(cmd)
        #combine cost and grad:
        cmd = 'paste cost.txt grad.txt > '+ file_name +' ; rm cost.txt grad.txt gradorig.txt'
        os.system(cmd)

        alist = open(file_name).read().split()
        iters  = numpy.asarray(alist[0::3])  
        cost = numpy.asarray(alist[1::3]).astype(np.float)
        grad = numpy.asarray(alist[2::3]).astype(np.float)
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
    plt.savefig('%s_%s.png'%(VAR[:4],expname[:-4]),dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
