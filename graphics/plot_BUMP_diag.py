import os 
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch

def readdata():
    '''
    BUMP log file name format: 
              mpas_*.info.0000
        mpas_3denvar_bump.info.0000  
        mpas_3dvar_bumpcov.info.0000  
        mpas_parametersbump_cov.info.0000

    '''
    bumplogfiles = []
    for files in os.listdir('../'):
        if fnmatch.fnmatch(files, 'mpas_*.info.0000'):
            bumplogfiles.append('../'+files)
    print(bumplogfiles)
    for bumplogfile in bumplogfiles:
        cmd = "sed -n -e '/--- Compute covariance/,/--- Compute correlation/p' \
              " + bumplogfile + " > blocks.txt"
        os.system(cmd)

        cmd = "sed -e '/--- Compute covariance/d' \
           -e '/Ensemble 1/d'        \
           -e '/Block 01_01_01_01/d' \
           -e '/Block 02_02_01_01/d' \
           -e '/Block 03_03_01_01/d' \
           -e '/Block 04_04_01_01/d' \
           -e '/Block 05_05_01_01/d' \
           -e '/----------------------/d' \
           -e '/--- Compute correlation/d' blocks.txt > tmp.txt"
        os.system(cmd)
        file_name = 'covariance_'+bumplogfile[3:][:-10]+'.txt'
        cmd = "awk '{print $2, $8}' tmp.txt > " + file_name + " ; rm blocks.txt tmp.txt"
        os.system(cmd)

        if (os.stat(file_name).st_size == 0):
            print(file_name+ " is empty.")
        else:
            lev = []
            cov = []
            with open(file_name,'r') as myfile:
 
                for line in myfile: #range(0, 54):  
                    all1 = numpy.asarray([a for a in filter(None, line.split(' '))])
                    lev.append(all1[0])
                    cov.append(all1[1])
                nlev = len(lev)/5
                print('Vertical Levels=',nlev) 
            plot(cov[0:nlev-1],lev[0:nlev-1],'Variance', 'T','($K^2$)',bumplogfile[3:][:-10])
            plot(cov[nlev:2*nlev-1],lev[nlev:2*nlev-1],'Variance', 'P', '($Pa^2$)',bumplogfile[3:][:-10])
            plot(cov[2*nlev:3*nlev-1],lev[2*nlev:3*nlev-1],'Variance', 'Q','($kg^2/kg^2$)',bumplogfile[3:][:-10])
            plot(cov[3*nlev:4*nlev-1],lev[3*nlev:4*nlev-1],'Variance', 'U','($m^2/s^2$)',bumplogfile[3:][:-10])
            plot(cov[4*nlev:5*nlev-1],lev[4*nlev:5*nlev-1],'Variance', 'V','($m^2/s^2$)',bumplogfile[3:][:-10])

def plot(var1,var2,var3,var4,var5,var6): 
    plt.plot(var1,var2,'g-*',markersize=5)
    plt.ylabel('Vertical Levels')
    plt.xlabel('%s %s'%(var3,var5))
    plt.legend('%s'%var4, loc='upper left')

    plt.grid(True)
    plt.savefig('cov_%s_%s.png'%(var4,var6),dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
