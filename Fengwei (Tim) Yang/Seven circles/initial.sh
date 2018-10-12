#$ -S /bin/bash
#$ -cwd
#$ -q mps.q
#$ -pe openmpi 1

##bash bash.sh

time ./PFBM.intel | tee output.txt

#time mpirun -np ${NSLOTS} ./PFBM.intel 1 4 | tee output.txt
