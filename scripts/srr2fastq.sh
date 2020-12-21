#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --output=%u.%x-%A.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --get-user-env
#SBATCH --mail-user=danny.bergeron@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=[0-20]

namelist=('SRR5637679.1'
'SRR5637680.1'
'SRR5637688.1'
'SRR5637689.1'
'SRR5637691.1'
'SRR6375597.1'
'SRR6375600.1'
'SRR6375596.1'
'SRR6375599.1'
'SRR6375595.1'
'SRR6375598.1'
'SRR6375593.1'
'SRR6375594.1'
'SRR6375592.1'
'SRR6375601.1'
'SRR5637466.1'
'SRR5637661.1'
'SRR5637467.1'
'SRR5637674.1'
'SRR5637470.1'
'SRR5637675.1'
)

name=${namelist[$SLURM_ARRAY_TASK_ID]}
fasterq-dump --skip-technical --split-files $name
