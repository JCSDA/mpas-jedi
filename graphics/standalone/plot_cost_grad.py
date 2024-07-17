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

    dalogfile = 'jedi.log'
    file_name = 'costgrad.txt'

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

    cmd = 'grep "'+VAR1+': Jb" '+ dalogfile +' \
        | grep -o -P "=\K.+"  > costJb.txt'
    os.system(cmd)

    cmd = 'grep "'+VAR1+': JoJc" '+ dalogfile +' \
        | grep -o -P "=\K.+"  > costJoJc.txt'
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
    cmd = 'paste cost.txt grad.txt costJb.txt costJoJc.txt > '+ file_name +' ; rm cost.txt grad.txt gradorig.txt costJb.txt costJoJc.txt'
    os.system(cmd)

    alist = open(file_name).read().split()
    iters  = numpy.asarray(alist[0::5])  
    cost = numpy.asarray(alist[1::5]).astype(np.float32)
    grad = numpy.asarray(alist[2::5]).astype(np.float32)
    costJb = numpy.asarray(alist[3::5]).astype(np.float32)
    costJoJc = numpy.asarray(alist[4::5]).astype(np.float32)

    forx=integers(1,len(iters)) 

    # plot(forx,cost,iters,VAR1,'cost')
    # plot(forx,grad,iters,VAR2,'grad')

    plot_in_one_figure(forx, cost, costJb, costJoJc, grad, iters)

def plot_in_one_figure(forx, cost, costJb, costJoJc, grad, iters):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Iterations',fontsize=12)
    ax1.set_xticks(forx[::10])
    ax1.set_xticklabels(iters.astype(np.int32)[::10]) #,rotation=45)
    ax1.set_ylabel('cost function',fontsize=12)
    plt.plot(forx,cost,'r.-', label='J')
    plt.plot(forx,costJb,'b.-', label='Jb')
    plt.plot(forx,costJoJc,'k.-', label='Jo')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,4))
    ax1.legend()

    ax2 = ax1.twinx()
    ax2.set_ylabel('gradient norm reduction',fontsize=16)
    ax2.plot(forx, grad, 'r.:', label='grad')
    #ax2.legend()

    plt.grid(True)
    plt.savefig('costgrad.png',dpi=200,bbox_inches='tight')
    plt.close()

def plot(forx,value,iters,VAR,name): 
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Iterations',fontsize=16)
    ax1.set_xticks(forx[::2])
    ax1.set_xticklabels(iters.astype(np.int32)[::2]) #,rotation=45)
    ax1.set_ylabel('%s'%VAR,fontsize=16)
    plt.plot(forx,value,'r.-')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,4))
    plt.grid(True)
    plt.savefig('%s.png'%(name),dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
