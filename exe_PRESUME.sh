 #!/bin/bash
#$ -S /bin/bash
#$ -N PRESUME
#$ -cwd

bash intermediate/shell/esu_${SGE_TASK_ID}.sh