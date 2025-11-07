# cd /lcrc/group/ATLAS/users/rwang/Argonne_computing/PPS-CCE/darshan/darshan_test/athena/
. /lcrc/project/ATLAS-HEP-group/rwang/Argonne_computing/PPS-CCE/darshan/darshan_test/athena/setup.sh
setupATLAS -c x86_64-el9-gcc13-opt  -m /lcrc:/lcrc --nohome -r "asetup 
Athena,main,latest;lsetup darshan;. venv_el9/bin/activate;python -m jupyterlab 
--no-browser --port 8080"
