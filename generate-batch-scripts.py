import pandas as pd
import argparse
import yaml

# Match tumour and normal sample accession codes based on sample name

class WESSample:

    sample_name = "submitted_subject_id"
    sample_type = "Sample Name"

    def __init__(self, name: str, filename: str, kw_labels: tuple):
        self.name = name
        self.runinfo_table = pd.read_csv(filename, sep=',', header=0) # TODO handle exception if file does not exist
        self.runs_per_patient = {}
        self.all_runs = self.runinfo_table["Run"].tolist()
        self.kw_to_label = {
            kw_labels[0] : "normal",
            kw_labels[1] : "tumour"
        }

    def preprocess_runinfo(self):
        # write/replace labels, assigning either "normal" or "tumour" to each run
        for pattern, standard_label in self.kw_to_label.items():
            # https://docs.python.org/3/library/re.html
            self.runinfo_table.replace({WESSample.sample_type: fr'.*{pattern}.*'}, {WESSample.sample_type: standard_label},  regex=True, inplace=True)

    def groupby_sample(self):
        # get all rows corresponding to a sample, create dictionary keys: sample and values: normal/tumour runs
        for sample in self.runinfo_table[WESSample.sample_name].unique():
            runs_per_type = dict.fromkeys(self.kw_to_label.values())
            for sample_type in self.kw_to_label.values():

                match_sample = self.runinfo_table[WESSample.sample_name] == sample
                match_sample_type = self.runinfo_table[WESSample.sample_type] == sample_type
                
                runs = self.runinfo_table.loc[match_sample & match_sample_type]["Run"].tolist()
                runs_per_type[sample_type] = runs
        
            self.runs_per_patient[sample] = runs_per_type

    @staticmethod
    def write_to_file(cmds: str, write_filename: str):
        with open (f"{output_directory}/{write_filename}", 'w') as fh:
            fh.write(cmds)

    # TODO change to switch/match case so can modify run to run.sra or prefix{run}.sra easily
    # maybe change parameter sra_suffix to accept both prefix and suffix
    def batch_reads_retrieval(self, scripts_directory: str, write_directory: str, sra_suffix = "", cmd="#!/bin/bash\n"):
        # TODO start from sra file, for restricted datasets
        cmd += "module load sratoolkit\n"
        for run in self.all_runs:
            run = run if sra_suffix == "" else f"{write_directory}/{run}{sra_suffix}.sra"
            cmd += (f"bash {scripts_directory}/fasterq-dump.sh {run} {write_directory}\n")
        WESSample.write_to_file(cmd.rstrip(), f"{name}-get-reads.sh")

    def batch_reads_quality_control(self, scripts_directory: str, raw_fastq_directory: str, write_directory: str, sra_suffix = "", cmd="#!/bin/bash\n"):
        cmd += "module load fastqc trimgalore\n"
        for run in self.all_runs:
            run = run if sra_suffix == "" else f"{run}{sra_suffix}" # HERE
            cmd += (f"bash {scripts_directory}/fastq-qc.sh {raw_fastq_directory}/{run}_1.fastq "
                f"{raw_fastq_directory}/{run}_1.fastq "
                f"{raw_fastq_directory}/{run}_2.fastq "
                f"{write_directory}\n")
        WESSample.write_to_file(cmd.rstrip(), f"{name}-quality-control-reads.sh")

    def batch_reads_mapping(self, scripts_directory: str, trimmed_fastq_directory: str, write_directory: str, sra_suffix = "", cmd="#!/bin/bash\n"):
        cmd += "module load bwa samtools\n"
        for run in self.all_runs:
            run = run if sra_suffix == "" else f"{run}{sra_suffix}" # HERE
            cmd += (f"{scripts_directory}/fastq-read-alignment.sh {trimmed_fastq_directory}/{run}_1_val_1.fq "
                f"{trimmed_fastq_directory}/{run}_2_val_2.fq {write_directory} {run}.bam\n")
        WESSample.write_to_file(cmd.rstrip(), f"{name}-mapping-reads.sh")

    def batch_bams_preprocessing(self, scripts_directory: str, mapped_reads_directory: str, write_directory: str, intervals_filename: str, sra_suffix = "", cmd="#!/bin/bash\n"):
        cmd += "module load GATK samtools\n"
        for run in self.all_runs:
            run = run if sra_suffix == "" else f"{run}{sra_suffix}" # HERE
            cmd += (f"bash {scripts_directory}/bam-process.sh {mapped_reads_directory}/{run}.bam {run}_rg "
                f"{write_directory} {run} {intervals_filename}\n")
        WESSample.write_to_file(cmd.rstrip(), f"{name}-preprocessing-bams.sh")

    def swarm_vcfs_mutation_calling(self, processed_bams_directory: str, write_directory: str, sra_suffix = "", cmd=""):
        _num = list(range(1,22+1))+['X','Y','M']
        chr = ["chr"+str(n) for n in _num]

        for sample, runs in self.runs_per_patient.items():
            # TODO multiple sites, this version only supports one paired samples
            normal_sample = runs["normal"][0]
            tumour_sample = runs["tumour"][0]
            sample_name = sample.replace(" ","")

            # HERE 
            normal_sample = normal_sample if sra_suffix == "" else f"{normal_sample}{sra_suffix}"
            tumour_sample = tumour_sample if sra_suffix == "" else f"{tumour_sample}{sra_suffix}"

            for chr_i in chr:
                cmd += ("module load GATK; gatk Mutect2 "
                    "-R /data/wuchh/somatic-mutation-calls/genomes/hg38/hg38.fa "
                    f"-L {chr_i} "
                    f"-I {processed_bams_directory}/{normal_sample}.recal.bam "
                    f"-I {processed_bams_directory}/{tumour_sample}.recal.bam "
                    f"-normal {normal_sample}_rg "
                    f"--germline-resource /data/wuchh/somatic-mutation-calls/databases/af-only-gnomad.hg38.vcf "
                    f"--panel-of-normals /data/wuchh/somatic-mutation-calls/databases/1000g_pon.hg38.vcf "
                    "--native-pair-hmm-threads $SLURM_CPUS_PER_TASK "
                    f"-O {write_directory}/{sample_name}_{chr_i}.vcf\n")
            WESSample.write_to_file(cmd.rstrip(), f"{name}-{sample_name}-mutect2.swarm")
            
    def batch_vcfs_annotations(self, processed_bams_directory: str, vcf_directory: str, write_directory: str, sra_suffix = "", cmd=""):
        _num = list(range(1,22+1))+['X','Y','M']

        # TODO fix this implementation  
        for sample, runs in self.runs_per_patient.items():
            normal_sample = runs["normal"][0] 
            tumour_sample = runs["tumour"][0]
            sample_name = sample.replace(" ","")

            # HERE 
            normal_sample = normal_sample if sra_suffix == "" else f"{normal_sample}{sra_suffix}"
            tumour_sample = tumour_sample if sra_suffix == "" else f"{tumour_sample}{sra_suffix}"

            vcf_per_chr = " ".join([f"I={vcf_directory}/{sample_name}_chr{str(n)}.vcf" for n in _num])
            cmd += f"java -Xmx16g -jar $PICARDJAR MergeVcfs {vcf_per_chr} O={write_directory}/merged_{sample_name}.vcf\n"

            vcf_stats_per_chr = " ".join([f"-stats {vcf_directory}/{sample_name}_chr{str(n)}.stats" for n in _num])
            cmd += f"gatk MergeMutectStats {vcf_per_chr} {vcf_stats_per_chr} -O {write_directory}/merged_{sample_name}.vcf.stats\n"

            cmd += (f"gatk GetPileupSummaries -I {processed_bams_directory}/{normal_sample}.recal.bam "
                f"-V /data/wuchh/somatic-mutation-calls/databases/af-only-gnomad.hg38.vcf "
                f"-L /data/wuchh/somatic-mutation-calls/databases/wgs_calling_regions.hg38.interval_list " # figure out which interval to use
                f"-O {vcf_directory}/normal_{sample_name}_pileup.table\n") 

            cmd += (f"gatk GetPileupSummaries -I {processed_bams_directory}/{tumour_sample}.recal.bam "
                f"-V /data/wuchh/somatic-mutation-calls/databases/af-only-gnomad.hg38.vcf "
                f"-L /data/wuchh/somatic-mutation-calls/databases/wgs_calling_regions.hg38.interval_list " # figure out which interval to use
                f"-O {vcf_directory}/tumour_{sample_name}_pileup.table\n") 

            cmd += (f"gatk CalculateContamination -I {vcf_directory}/tumour_{sample_name}_pileup.table "
                f"-matched {vcf_directory}/normal_{sample_name}_pileup.table "
                f"-O {vcf_directory}/{sample_name}_contamination.table\n")

            cmd += ("gatk FilterMutectCalls -R /data/wuchh/somatic-mutation-calls/genomes/hg38/hg38.fa "
                f"-V merged_{sample_name}.vcf "
                f"--contamination-table {vcf_directory}/{sample_name}_contamination.table "
                f"-O {vcf_directory}/filtered_{sample_name}.vcf\n")

            cmd += ("java -Xmx16g -jar $SNPEFF_JAR -noStats -hgvs -lof hg38 "
                f"{vcf_directory}/filtered_{sample_name}.vcf > {vcf_directory}/annotated_{sample_name}.vcf")

            WESSample.write_to_file(cmd.rstrip(), f"{name}-{sample_name}-filter-annotate.sh")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--name')
    parser.add_argument('-I', '--infile')
    parser.add_argument('-O', '--outdir')
    parser.add_argument('--normal')
    parser.add_argument('--tumour')
    parser.add_argument('--config')
    parser.add_argument('-L', '--intervals', nargs='?', default='/data/wuchh/somatic-mutation-calls/databases/wgs_calling_regions.hg38.interval_list')
    args = parser.parse_args()

    name = args.name
    metadata_filename = args.infile
    output_directory = args.outdir
    kw_normal, kw_tumour = args.normal, args.tumour
    config_filename = args.config
    intervals_filename = args.intervals

    with open(config_filename, mode="rb") as file:
        file_paths = yaml.safe_load(file)

    scripts = file_paths["directories"]["precall_scripts"]
    raw_reads = file_paths["directories"]["raw_reads"]
    trimmed_reads = file_paths["directories"]["trimmed_reads"]
    aligned_reads = file_paths["directories"]["mapped_reads"]
    processed_bams = file_paths["directories"]["bam_files"]
    mutation_vcf = file_paths["directories"]["variant_call_files"]
    
    example = WESSample(name, metadata_filename, (kw_normal, kw_tumour))
    example.preprocess_runinfo()
    example.groupby_sample()
    example.batch_reads_retrieval(scripts, raw_reads, sra_suffix="_dbGaP-25281")
    example.batch_reads_quality_control(scripts, raw_reads, trimmed_reads, sra_suffix="_dbGaP-25281")
    example.batch_reads_mapping(scripts, trimmed_reads, aligned_reads, sra_suffix="_dbGaP-25281")
    example.batch_bams_preprocessing(scripts, aligned_reads, processed_bams, intervals_filename, sra_suffix="_dbGaP-25281")
    example.swarm_vcfs_mutation_calling(processed_bams, mutation_vcf, sra_suffix="_dbGaP-25281")
    example.batch_vcfs_annotations(processed_bams, mutation_vcf, mutation_vcf, sra_suffix="_dbGaP-25281")

    # TODO organize output files
