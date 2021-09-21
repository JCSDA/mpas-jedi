#!/bin/csh
echo "test cylc setup" 
module purge
module load cylc
module load graphviz

setenv TOP_DIR  /gpfs/fs1/p/mmm/parc/${USER}/test/pandac/work
setenv QUEUE    regular
setenv ACC_KEY  NMMM0035
setenv currdir  `basename "$PWD"`
echo 'currdir=' $currdir

rm -fr ${HOME}/cylc-run/$currdir
mkdir /glade/scratch/${USER}/pandac/cylc-run_$currdir


cat >! suite.rc << EOF
#!Jinja2
{% set PARAMS = range(0,11) %}
title = "cycling jobs"
[cylc]
    UTC mode = False
    [[environment]]
[scheduling]
    max active cycle points = 200
    initial cycle point = 20180415T00
    final cycle point   = 20180417T00
    [[dependencies]]
        # Initial cycle point
        [[[R1]]]    # Run once, at the initial point.
            graph = "DA => MPASFC1"
        [[[PT6H]]]
            graph = "MPASFC1[-PT6H] => DA => MPASFC1"
        [[[P1D]]]
            graph = "DA => MPASFC2"
        [[[P1D]]]
            graph =   """
            {% for I in PARAMS%}
               MPASFC2 => OMF{{I}}
            {% endfor %}
                      """
        [[[PT6H]]]
            graph = "DA => DADIAG"
        [[[P1D]]]
            graph = "MPASFC2 => FCDIAG"
[runtime]
    [[root]] # suite defaults
        pre-script = "cd  ${TOP_DIR}/$currdir/"
    [[DA]]
        script = ${TOP_DIR}/$currdir/da.csh \${CYLC_TASK_CYCLE_POINT}
        [[[job submission]]]
                method = pbs
                execution time limit = PT29M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=4:ncpus=9:mpiprocs=9
                -q = ${QUEUE}
                -A = ${ACC_KEY}
    [[MPASFC1]]
        script = ${TOP_DIR}/$currdir/fc1.csh
        [[[job submission]]]
                method = pbs
                execution time limit = PT9M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=4:ncpus=32:mpiprocs=32
                -q = ${QUEUE}
                -A = ${ACC_KEY}
    [[MPASFC2]]
        script = ${TOP_DIR}/$currdir/fc2.csh
        [[[job submission]]]
                method = pbs
                execution time limit = PT1H10M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=4:ncpus=32:mpiprocs=32
                -q = ${QUEUE}
                -A = ${ACC_KEY}
{% for PARAM  in PARAMS %}
    [[OMF{{PARAM}}]]
        script = ${TOP_DIR}/$currdir/omf.csh {{PARAM}}
        [[[job submission]]]
                method = pbs
                execution time limit = PT3M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=4:ncpus=9:mpiprocs=9
                -q = ${QUEUE}
                -A = ${ACC_KEY}
{% endfor %}
    [[DADIAG]]
        script = ${TOP_DIR}/$currdir/dadiag.csh
        [[[job submission]]]
                method = pbs
                execution time limit = PT40M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=1:ncpus=32:mpiprocs=32
                -q = ${QUEUE}
                -A = ${ACC_KEY}
    [[FCDIAG]]
        script = ${TOP_DIR}/$currdir/fcdiag.csh
        [[[job submission]]]
                method = pbs
                execution time limit = PT29M
        [[[directives]]]
                -j = oe
                -S = /bin/csh
                -l = select=1:ncpus=1:mpiprocs=1
                -q = ${QUEUE}
                -A = ${ACC_KEY}
[visualization]
    initial cycle point = 20180415T00
    final cycle point   = 20180417T00
    number of cycle points = 20
    default node attributes = "style=filled", "fillcolor=grey"
EOF

cylc register $currdir  /gpfs/fs1/p/mmm/parc/${USER}/test/pandac/work/$currdir
cylc validate $currdir
cylc run $currdir 

exit
