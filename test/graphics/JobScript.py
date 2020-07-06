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

    Each HPC job submission system (e.g., PBSPro on Cheyenne)
    will have its own derived class that defines
    the job header and submission command
    generic config elements:
        required config parameter(s):
        basescript (required) - either a list of str's containing individual lines of the script
                                or a str giving the location of the script

        optional config parameter(s):
        env - linux environment of the script (e.g., csh, bash, sh, tcsh)
        name - job name
        nppernode - processors per node
        nnode - number of nodes
        walltime - walltime
        olog - output log name
        elog - error log name
    '''
    def __init__(self, conf):
        ## job descriptors
        self.name = conf.get('name','PyJobScript')
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
        self.script = self.name+'.job.'+self.env

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
        self.account = conf.get('account','NMMM0015')
        self.queue = conf.get('queue','regular')
        assert self.queue in self.qavail, ("ERROR: PBSProCheyenne requires queue to be any of ",self.qavail)
        self.memory = conf.get('memory',109)
        assert self.memory in self.mavail, ("ERROR: PBSProCheyenne requires memory (in GB) to be any  of", self.mavail)
        assert self.nppernode <= self.maxnppernode, ("ERROR: PBSProCheyenne requires nppernode <= ", self.maxnppernode)

        self.header = [
            '#PBS -N '+self.name,
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


JobScriptDict = {
    ## cheyenne
    # login nodes
    'cheyenne': PBSProCheyenne,
    # cron jobs
    'chadmin': PBSProCheyenne,
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
