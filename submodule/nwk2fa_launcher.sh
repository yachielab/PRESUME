 #!/bin/bash
#$ -S /bin/bash
#$ -N PRESUME
#$ -cwd

SCRIPT_PATH=$1

bash ${SCRIPT_PATH}/downstream_${SGE_TASK_ID}.sh