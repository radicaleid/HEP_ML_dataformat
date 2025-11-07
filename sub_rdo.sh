#!/bin/bash
export START_TIME=$(date +'%Y-%m-%d %H:%M:%S')

rel=main,latest #main--dev3LCG,r02-11 #25.0.25 #25.0.11 #23.0.53 # 24.0.41 #24.0.23 #'23.0.53' #Athena,23.0.53 # for workflow
export release=Athena,$rel

########### RDO MRG #####################################
proc='ttbar'
athena='overlay-rec'
sharedWriter=false
parallelCompression=false

stage='RDODump'  #'RDOMerge'
workdir=workdir/test/RDODump/${rel//,/-}/$proc

for format in 'RNTuple' #'TTree' #'RNTuple'
do
    for nevt in 50 #100 500 1000 1200 2400 3600
    do
        for nproc in  2
        do
            # setupATLAS -c el9 -m /lcrc:/lcrc --nohome -r " export release=${release} && /srv/run_rec.sh $athena $nevt $nproc $sharedWriter $proc $format $parallelCompression Reconstruction/RNTuple/${rel//,/_}/$proc"
            # args='setupATLAS -c el9 -m /lcrc:/lcrc -r " export release=${release} && /srv/run_FastChain.sh $stage $nevt $nproc $lessCompression $proc $format $parallelCompression workdir/FastChain/${rel//,/-}/$proc" '
            # sbatch -c $nproc -t 00:29:00 ./submit_singularity.sh $args # sbatch job
            # echo "$args" | bash
            # setupATLAS -c el9 -m /lcrc:/lcrc -r " export release=${release} &&  && /srv/run_FastChain.sh $stage $nevt $nproc $lessCompression $proc $format $parallelCompression workdir/FastChain/test/Track_Overlay"
           /cvmfs/atlas.cern.ch/repo/containers/sw/apptainer/x86_64-el7/1.2.2/bin/apptainer exec -B /gpfs:/gpfs -B /lcrc:/lcrc -B /cvmfs:/cvmfs -B $PWD:/srv /cvmfs/atlas.cern.ch/repo/containers/fs/singularity/x86_64-almalinux9 /srv/run_rdo.sh $stage $nevt $nproc $sharedWriter $proc $format $parallelCompression $workdir

        done
    done
done