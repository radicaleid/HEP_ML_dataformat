#!/bin/bash
localdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

workflow=$1
nevents_per_proc=$2 #5
nproc=$3 #32
sharedWriter=$4 # TrackOverlay: lessCompression
proc=$5
format=$6 #'MCOverlay' #'TrackOverlay'
parallelCompression=$7
workdir=$localdir/$8

echo $@
echo $workflow $nevents_per_proc $nproc $sharedWriter $proc $format $parallelCompression $workdir

darshan_setup_file=/lcrc/group/ATLAS/users/rwang/Argonne_computing/PPS-CCE/darshan/darshan_test/athena/configs/FastChain_lcrc_cvmfs.sh
darshan_conf_file=/lcrc/group/ATLAS/users/rwang/Argonne_computing/PPS-CCE/darshan/darshan_test/athena/configs/FastChain_wf1_sim_env_TACC.conf
echo [$SECONDS]create working dir $workdir
mkdir -p $workdir
echo [$SECONDS]copy darshan setup to $workdir
cp -v $darshan_setup_file $workdir/setup.sh
chmod 777 $workdir/setup.sh
echo [$SECONDS] copy $darshan_conf_file to $workdir
cp -v $darshan_conf_file $workdir/darshan.conf
chmod 777 $workdir/darshan.conf
echo -e 'asetup $release' >> $workdir/setup.sh

case ${proc} in
"ttbar")
    RDO_File='/lcrc/group/ATLAS/users/rwang/Argonne_computing/data/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8481_s4149_r15238/*'
    if [ $format = 'MCOverlay' ]
    then 
        # RDO_BKG_File=$localdir/configs/RDO.33837529._000247.pool.root.1
        RDO_BKG_File=/lcrc/group/ATLAS/users/rwang/Argonne_computing/data/valid2.900149.PG_single_nu_Pt50.merge.RDO.e8514_e8528_s4332_d1960_d1943_tid40126453_00/* # MC23e
        
    else
        # RDO_BKG_File=$localdir/configs/RDO_BKG_10k_MC23c.pool_0.root
        # RDO_BKG_File=$localdir/configs/data_prep/RDO_BKG_MC23/mc23_13p6TeV.900149.PG_single_nu_Pt50.digit.RDO.e8514_e8528_s4153_d1879/RDO.*.pool.root.1
        RDO_BKG_File=/lcrc/group/ATLAS/users/rwang/Argonne_computing/data/valid1.900149.PG_single_nu_Pt50.recon.RDO.e8514_e8528_s4332_d1960_d1943_r16231_tid42346764_00/*
#$localdir/configs/RDO_BKG_2k.pool_0.root
    fi
esac
RDO_MRGFile="RDO.pool.root"
ESD_File="ESD.pool.root"
AOD_File="AOD.pool.root"
NTUP_File="NTUP.pool.root"

# export DARSHAN_DEFAULT_NPROCS=7
events=$(($nproc * $nevents_per_proc))
subfolder=$(date +'%Y/%m/%d')
logfolder=$DARSHAN_LOG_DIR/${subfolder//"/0"/"/"}

case ${workflow} in
"RDOMerge")
    #--- athenaMP ---
    echo [$SECONDS]$format
    workdir=$workdir/${workflow}_${proc}_${format}_${nproc}_${nevents_per_proc}_${sharedWriter}_${parallelCompression}
    cmd=''
    preExec_cmd=''
    postExec_cmd=''
    if [[ $format =~ RNTuple ]]; then
        preExec_cmd+="flags.Output.StorageTechnology.EventData=\"ROOTRNTUPLE\";"
    fi
    echo [$SECONDS]"working in $workdir"
    (rm -rf $workdir && mkdir -p $workdir)
    ! (cd $workdir && export ATHENA_CORE_NUMBER=${nproc} && . ../setup.sh \
    && RDOMerge_tf.py \
    --maxEvents ${events} --skipEvents 0 \
    --multiprocess --athenaMPMergeTargetSize "RDO*:0" --outputFileValidation False \
    --inputRDOFile ${RDO_File} --outputRDO_MRGFile ${RDO_MRGFile} \
    --postInclude "default:PyJobTransforms.UseFrontier" \
    --preExec ${preExec_cmd} \
    --athenaopts "  --preloadlib=$ATLASMKLLIBDIR_PRELOAD/libintlc.so.5:$ATLASMKLLIBDIR_PRELOAD/libimf.so:$DARSHAN_LD_PRELOAD" \
    --imf False ${cmd} \
    && echo "Job end:" $(date +'%Y-%m-%d %H:%M:%S')$'\n' \
    # && . $localdir/venv_el9/bin/activate && python3 -m darshan summary --enable_dxt_heatmap 20*/*/*/*.darshan \
    )  2>&1 |tee $workdir.log
;;
"RDODump")
    #--- athenaMP ---
    echo [$SECONDS]$format
    workdir=$workdir/${workflow}_${proc}_${format}_${nproc}_${nevents_per_proc}_${sharedWriter}_${parallelCompression}
    cmd=''
    preExec_cmd=''
    postExec_cmd=''
    if [[ $format =~ RNTuple ]]; then
        preExec_cmd+="--preExec 'flags.Output.StorageTechnology.EventData=\"ROOTRNTUPLE\";'"
        RDO_File=${workdir//RDODump/RDOMerge}/RDO*root*
        echo [$SECONDS]input: ${RDO_File}
    fi
    echo [$SECONDS]"working in $workdir"
    (rm -rf $workdir && mkdir -p $workdir)
    ! (cd $workdir  && . ../setup.sh \
    && ATHENA_CORE_NUMBER=${nproc} Reco_tf.py \
    --CA 'all:True' --autoConfiguration 'everything' \
    --conditionsTag 'all:OFLCOND-MC21-SDR-RUN4-02' \
    --geometryVersion 'all:ATLAS-P2-RUN4-03-00-00' \
    --multiprocess --athenaMPMergeTargetSize "AOD*:0" --outputFileValidation False \
    --steering 'doRAWtoALL' \
    --digiSteeringConf 'StandardInTimeOnlyTruth' \
    --postInclude 'all:PyJobTransforms.UseFrontier,InDetConfig.SiSpacePointFormationConfig.InDetToXAODSpacePointConversionCfg' \
    --preInclude 'all:Campaigns.PhaseIIPileUp200' 'InDetConfig.ConfigurationHelpers.OnlyTrackingPreInclude'\
    --postExec 'from InDetGNNTracking.InDetGNNTrackingConfig import DumpObjectsCfg; cfg.merge(DumpObjectsCfg(flags));cfg.getService(\"AthMpEvtLoopMgr\").ExecAtPreFork=[\"AthCondSeq\"]; cfg.getService(\"AthenaPoolCnvSvc\").PoolAttributes += [ \"DatabaseName = '*'; RNTUPLE_WRITER_METRICS_ENABLED = '1'; \"];'\
    ${preExec_cmd} \
    --maxEvents ${events} --skipEvents 0 \
    --inputRDOFile ${RDO_File} \
    --outputESDFile ${ESD_File} \
    --outputAODFile 'test.aod.gnnreader.debug.root' \
    --athenaopts "  --preloadlib=$ATLASMKLLIBDIR_PRELOAD/libintlc.so.5:$ATLASMKLLIBDIR_PRELOAD/libimf.so:$DARSHAN_LD_PRELOAD" \
    --imf False ${cmd} \
    && echo "Job end:" $(date +'%Y-%m-%d %H:%M:%S')$'\n' \
    # && . $localdir/venv_el9/bin/activate && python3 -m darshan summary --enable_dxt_heatmap 20*/*/*/*.darshan \
    )  2>&1 |tee $workdir.log
esac

# --multithreaded=True 
    # --sharedWriter ${sharedWriter} --parallelCompression ${parallelCompression} \
########## write root ntyple for training ##########
    # --postExec 'from InDetGNNTracking.InDetGNNTrackingConfig import DumpObjectsCfg; cfg.merge(DumpObjectsCfg(flags))' \
########## debug message #############
# --postExec 'cfg.getService("MessageSvc").defaultLimit=999999;cfg.getService("MessageSvc").OutputLevel=0;