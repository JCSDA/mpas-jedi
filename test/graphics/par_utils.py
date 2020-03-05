import getopt
import os

#================================
# parallel processing definitions
#================================

def proc_print(myproc,msg):
    # addendum to normal print()
    # add prefix of processor rank (myproc) when myproc > 0
    prefix = ""
    if myproc > 0:
        prefix = "p="+"{:d}".format(myproc)+": "
    print(prefix+msg)

def print_par_args_usage(thisPyScript):
    print ('Either 0 or 2 arguments are required.')
    print ('usage: python '+thisPyScript+'')
    print ('   or  python '+thisPyScript+' -n <NP> -i <proc>')
    print ('   or  python '+thisPyScript+' --nprocs <NP> --iproc <proc>')
    print ('')
    print ('  where <NP> is the total number of processors')
    print ('  and <proc> is the processor rank, starting at 1')
    print ('')
    print('   E.g., using GNU parallel:')
    print ('')
    print('   parallel -j<NP> --plus "python '+thisPyScript+' -n {##} -i {#} >& diags{#}.log" ::: `seq <NP>`')
    os._exit(-1)

def par_args(argv):
    opts, args = getopt.getopt(argv[1:], 'n:i:', ['nprocs', 'iproc'])
    if len(opts) == 0:
        return [1,0]
    elif len(opts) == 2:
        nproc = int(opts[0][1])
        myproc = int(opts[1][1]) - 1
        if myproc < 0:
           print_par_args_usage(argv[0])
        return [nproc, myproc]
    else:
        print_par_args_usage(argv[0])

#================================
#================================

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
