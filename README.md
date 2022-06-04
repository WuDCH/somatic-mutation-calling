# somatic-mutation-calling
Pipeline to call somatic variants from matched tumour-normal pairs

**CPU and memory allocation for each step**
sbatch --time=24:00:00 --gres:lscratch:256 --mem=32g --cpus-per-task=16 <fasterq-dump.sh>
sbatch --time=24:00:00 --gres:lscratch:128 --mem=32g <trim.sh>
