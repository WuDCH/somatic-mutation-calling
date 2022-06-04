# somatic-mutation-calling
Pipeline to call somatic variants from matched tumour-normal pairs

**CPU and memory allocation for each step**

sbatch --time=24:00:00 --gres:lscratch:256 --mem=32g --cpus-per-task=16 <fasterq-dump.sh>

sbatch --time=24:00:00 --gres:lscratch:128 --mem=32g <trim.sh>

sbatch --time=48:00:00 --mem=64g --cpus-per-task=32 <alignment.sh>

sbatch --time=48:00:00 --gres=lscratch:500 --mem=64g --cpus-per-task=32 <alignment.sh>

python3 generate-batch-scripts.py --name VanAllen-2015 -I ~/Downloads/VanAllen_SRAruntable.csv -O /Users/dch/Downloads --normal normal --tumour tumor --config /Users/dch/Documents/somatic-mutation-calling/sample.yaml -L /data/wuchh/somatic-mutation-calls/databases/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed
