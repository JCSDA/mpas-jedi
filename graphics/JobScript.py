#!/bin/env python3

import os
from pathlib import Path
import subprocess

class JobScriptBase():
    '''
    Describes an HPC job script including
    + reading a configuration
    + generating the script
    + submitting the script

    Each HPC job submission system (e.g., PBSPro on Cheyenne/Casper/Derecho)
    will have its own derived class that defines
    the job header and submission command
    generic config elements:
        required config parameter(s):
        basescript (required) - either a list of str's containing individual lines of the script
                                or a str giving the location of the script

        optional config parameter(s):
        env - linux environment of the script (e.g., csh, bash, sh, tcsh)
        name - job name
        filename - job script file name
        nppernode - processors per node
        nnode - number of nodes
        walltime - walltime
        olog - output log name
        elog - error log name
    '''
    def __init__(self, conf):
        ## job descriptors
        self.jobname = conf.get('name', 'PyJobScript')
        self.filename = conf.get('filename', self.jobname)
        self.nppernode = conf.get('nppernode',1)
        self.nnode = conf.get('nnode',1)
        self.walltime = conf.get('walltime','01:00:00')
        self.olog = conf.get('olog','log.job.out')
        self.elog = conf.get('elog','log.job.err')

        ## submission descriptors
        self.env = conf.get('env','csh')
        self.basescript = conf['script']
        assert (isinstance(self.basescript,list) or isinstance(self.basescript,str)), \
            "JobScriptBase : basescript must either be a list or a string"
        self.jobpath = Path(conf.get('path','./'))
        self.script = self.filename+'.job.'+self.env

        self.command = './'
        self.header = []

    def create(self):
        ## initialize env
        joblines = ['#!/bin/'+self.env+'\n']

        ## concatenate header
        for line in self.header: joblines.append(line+'\n')

        ## concatenate script
        if isinstance(self.basescript,list):
            # assume self.basescript is a list of lines
            for line in self.basescript: joblines.append(line)

        elif isinstance(self.basescript,str):
            # assume self.basescript is a file
            bs = open(self.basescript,'r')
            for line in bs: joblines.append(line)
            bs.close()

        ## create the job path
        self.jobpath.mkdir(parents=True, exist_ok=True)

        ## write job script
        script = str(self.jobpath/self.script)
        if os.path.exists(script):
            os.remove(script)
        js = open(script,'w')
        js.writelines(joblines)
        js.close()
        os.system('chmod 744 '+script)

    def submit(self):
        ## submit job
        command = self.command+self.script
        CWD = os.getcwd()
        os.chdir(str(self.jobpath))
        print(command+" in "+os.getcwd())
        os.system(command)
        os.chdir(CWD)

class PBSProDerecho(JobScriptBase):
    '''
    PBSPro job script on Derecho
    unique config elements compared to base class:
        account - derecho account for charging
        queue   - name of job submission queue (see qavail)
        memory  - amount of memory requested per node (see mavail)

    NOTE: Derecho has a maximum of 128 processors available per node
    '''
    qavail = ['main']
    mavail = [235]
    maxnppernode = 128
    def __init__(self, conf):
        # Initialize derived config settings
        super().__init__(conf)

        # Initialize config settings that are specific to PBSProDerecho
        self.account = conf.get('account','NMMM0043')
        self.queue = conf.get('queue','main')
        assert self.queue in self.qavail, ("ERROR: PBSProDerecho requires queue to be any of ",self.qavail)
        self.memory = conf.get('memory',235)
        assert self.memory in self.mavail, ("ERROR: PBSProDerecho requires memory (in GB) to be any  of", self.mavail)
        assert self.nppernode <= self.maxnppernode, ("ERROR: PBSProDerecho requires nppernode <= ", self.maxnppernode)

        self.header = [
            '#PBS -N '+self.jobname,
            '#PBS -A '+self.account,
            '#PBS -q '+self.queue,
            '#PBS -l select='+str(self.nnode)+':ncpus='+str(self.nppernode)+':mpiprocs='+str(self.nppernode)+':mem='+str(self.memory)+'GB',
            '#PBS -l walltime='+self.walltime,
            '#PBS -m ae',
            '#PBS -k eod',
            '#PBS -o '+self.olog,
            '#PBS -e '+self.elog,
        ]

        self.command = 'qsub '


class PBSProCheyenne(JobScriptBase):
    '''
    PBSPro job script on Cheyenne
    unique config elements compared to base class:
        account - cheyenne account for charging
        queue   - name of job submission queue (see qavail)
        memory  - amount of memory requested per node (see mavail)

    NOTE: Cheyenne has a maximum of 36 processors available per node
    '''
    qavail = ['economy', 'regular', 'premium']
    mavail = [45, 109]
    maxnppernode = 36
    def __init__(self, conf):
        # Initialize derived config settings
        super().__init__(conf)

        # Initialize config settings that are specific to PBSProCheyenne
        self.account = conf.get('account','NMMM0043')
        self.queue = conf.get('queue','regular')
        assert self.queue in self.qavail, ("ERROR: PBSProCheyenne requires queue to be any of ",self.qavail)
        self.memory = conf.get('memory',109)
        assert self.memory in self.mavail, ("ERROR: PBSProCheyenne requires memory (in GB) to be any  of", self.mavail)
        assert self.nppernode <= self.maxnppernode, ("ERROR: PBSProCheyenne requires nppernode <= ", self.maxnppernode)

        self.header = [
            '#PBS -N '+self.jobname,
            '#PBS -A '+self.account,
            '#PBS -q '+self.queue,
            '#PBS -l select='+str(self.nnode)+':ncpus='+str(self.nppernode)+':mpiprocs='+str(self.nppernode)+':mem='+str(self.memory)+'GB',
            '#PBS -l walltime='+self.walltime,
            '#PBS -m ae',
            '#PBS -k eod',
            '#PBS -o '+self.olog,
            '#PBS -e '+self.elog,
        ]

        self.command = 'qsub '


class SLURMCasper(JobScriptBase):
    '''
    SLURM job script on Casper
    unique config elements compared to base class:
        account - casper account for charging
        partition - name of casper partition (see pavail)
        memory  - amount of memory requested per node (see maxmemory)

    NOTE: Casper has a maximum of 36 processors available per node
    '''
    pavail = ['dav']
    maxnppernode = 36
    maxmemory = 1100
    def __init__(self, conf):
        # Initialize derived config settings
        super().__init__(conf)

        # Initialize config settings that are specific to SLURMCasper
        self.account = conf.get('account','NMMM0015')
        self.partition = conf.get('partition','dav')
        assert self.partition in self.pavail, ("ERROR: SLURMCasper requires partition to be any of ",self.pavail)
        self.memory = conf.get('memory',300)
        assert self.memory <= self.maxmemory, ("ERROR: SLURMCasper requires memory (in GB) to be <= ", self.maxmemory)

        assert self.nppernode <= self.maxnppernode, ("ERROR: SLURMCasper requires nppernode <= ", self.maxnppernode)

        self.header = [
            '#SBATCH --job-name='+self.jobname,
            '#SBATCH --account='+self.account,
            '#SBATCH --ntasks='+str(self.nnode),
            '#SBATCH --cpus-per-task='+str(self.nppernode),
            '#SBATCH --mem='+str(self.memory)+'G',
            '#SBATCH --time='+self.walltime,
            '#SBATCH --partition='+self.partition,
            '#SBATCH --output='+self.olog,
        ]

        self.command = 'sbatch '


class PBSProCasper(JobScriptBase):
    '''
    PBSPro job script on Casper
    unique config elements compared to base class:
        account - casper account for charging
        queue   - name of job submission queue (see qavail)
        memory  - amount of memory requested per node (see mavail)

    NOTE: 96 of the 104 Casper compute nodes only have 36 processors. 8 GPU nodes have 128 processors.
    '''
    qavail = ['casper', 'gpudev'] # 'casper' queue is the default submission queue
    maxnppernode = 36
    maxmemory = 360
    def __init__(self, conf):
        # Initialize derived config settings
        super().__init__(conf)

        # Initialize config settings that are specific to PBSProCasper
        self.account = conf.get('account','NMMM0015')
        self.queue = conf.get('queue','casper')
        assert self.queue in self.qavail, ("ERROR: PBSProCasper requires queue to be any of ",self.qavail)
        # to use devices with gpus, use the 'casper' queue but specify ngpus in the select statement
        if self.queue == 'gpudev':
          self.queue = 'casper'
          gpus = ':ngpus=1'
        else:
          gpus = ""

        self.memory = conf.get('memory',109)
        assert self.memory <= self.maxmemory, ("ERROR: PBSProCasper requires memory (in GB) to be <= ", self.maxmemory)
        if self.nppernode > self.maxnppernode:
          print("WARNING: PBSProCasper requires nppernode <= ", self.maxnppernode)
          print(" changing nppernode from ", self.nppernode, " to ", self.maxnppernode)
          self.nppernode = self.maxnppernode
        assert self.nppernode <= self.maxnppernode, ("ERROR: PBSProCasper requires nppernode <= ", self.maxnppernode)

        self.header = [
            '#PBS -N '+self.jobname,
            '#PBS -A '+self.account,
            '#PBS -q '+self.queue,
            '#PBS -l select='+str(self.nnode)+':ncpus='+str(self.nppernode)+':mpiprocs='+str(self.nppernode)+':mem='+str(self.memory)+'GB'+gpus,
            '#PBS -l walltime='+self.walltime,
            '#PBS -m ae',
            '#PBS -k eod',
            '#PBS -o '+self.olog,
            '#PBS -e '+self.elog,
        ]

        self.command = 'qsub '


JobScriptDict = {
    ## derecho
    # computing nodes
    'derecho': PBSProDerecho,
    ## cheyenne
    # login nodes
    'cheyenne': PBSProCheyenne,
    # cron jobs
    'chadmin': PBSProCheyenne,
    ## casper
    'casper': PBSProCasper
}


def JobScriptFactory(conf):
    ## get system name
    conf['sysname'] = subprocess.run(['uname','-n'],
                      stdout=subprocess.PIPE).stdout.decode('utf-8')

    ## match system name with JobScriptDict or return base class object by default
    for key, jobclass in JobScriptDict.items():
        if key in conf['sysname']:
            return jobclass(conf)
    return JobScriptBase(conf)
