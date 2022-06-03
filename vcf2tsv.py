from numpy import int16, int8
import pandas as pd
import os
import re

class VCF_UTILS:

    def __init__(self, vcf_file_folder):
        self.vcf_directory = vcf_file_folder
        self.formatted_vcf_collector = []

    @staticmethod
    def find_starting_line(file, skiplines=0):
        # find start of line
        with open (file) as file_handle:
            for line in file_handle:
                if not line.startswith("##"):
                     return skiplines
                skiplines +=  1

        return ("LINE NOT FOUND") # throw exception

    def write_format(df, idx_format_col, format_labels, stype):
        _read_info = df.iloc[:,idx_format_col].str.split(pat=":", expand=True).iloc[:,:len(format_labels)]
        _read_info.columns = format_labels

        _read_info["BIALLELIC"] = _read_info["AD"].apply(lambda x: x.count(',') == 1)
        _read_info = _read_info[_read_info["BIALLELIC"] == True] # for now ignore multi-allelic

        _counts = _read_info["AD"].str.split(pat=",", expand=True).astype(int16)
        _counts.columns = [stype+"_REF_READS", stype+"_ALT_READS"]
        _counts[stype+"_TOTAL_READS"] = _read_info["DP"]

        return _counts

    @staticmethod
    def parse_vcf_file(vcf_file, pat_id):

        # find start of line
        _start = VCF_UTILS.find_starting_line(vcf_file)

        vcf_to_dataframe = pd.read_csv(vcf_file, skiprows=_start, sep='\t')

        base = [
            "#CHROM",
            "POS",
            "REF",
            "ALT"
        ]

        base_columns = vcf_to_dataframe.loc[:,base]

        _info = [
            "Allele",
            "Annotation",
            "Annotation_Impact",
            "Gene_Name",
            "Gene_ID",
            "Feature_Type",
            "Feature_ID",
            "Transcript_BioType",
            "Rank",
            "HGVS.c",
            "HGVS.p"
        ]
        print(vcf_to_dataframe)

        # in instances where there are multiple annotations for one variant, we consider the first
        # as being representative of the most deleterious change.
        # annotations are delimited by the keyword "ANN", collect first n entries after

        _variant_annotations = vcf_to_dataframe["INFO"].apply(lambda x: x[x.find("ANN"):])
        _variant_annotations = _variant_annotations.str.split(pat="|", n=len(_info), expand=True).iloc[:,:len(_info)]
        _variant_annotations.columns =_info
        _variant_annotations.drop(["Allele","Annotation_Impact","Gene_ID","Rank"], axis=1, inplace=True)

        _format = [
            "GT",
            "AD",
            "AF",
            "DP"
        ]

        # ALT REF

        idx_normal = -2
        idx_baseline = -1

        normal_read_columns = VCF_UTILS.write_format(vcf_to_dataframe, idx_normal, _format, "normal")
        baseline_read_columns = VCF_UTILS.write_format(vcf_to_dataframe, idx_baseline, _format, "baseline")
       
        idx_biallelic = list(normal_read_columns.index) 
        _variant_annotations = _variant_annotations.reindex(idx_biallelic)
        base_columns = base_columns.reindex(idx_biallelic)

        # add sample column
        samples_column = [pat_id] * len(idx_biallelic)


        _df_to_merge = [
            base_columns,
            _variant_annotations,
            normal_read_columns,
            baseline_read_columns,
        ]

        formatted_vcf = pd.concat(_df_to_merge, axis=1)
        formatted_vcf["Sample_ID"] = samples_column
        
        return formatted_vcf

    def merge_vcdf(self):
        for file in os.listdir(self.vcf_directory):
            p_id = (re.split('_',file)[-1].split(".")[0])
            f_vcf = VCF_UTILS.parse_vcf_file(f"{self.vcf_directory}/{file}", p_id)
            self.formatted_vcf_collector.append(f_vcf)

        concat_vcf = pd.concat(self.formatted_vcf_collector, axis=0)
        return concat_vcf

if __name__ == "__main__":
    
    a = VCF_UTILS("/Users/dch/Downloads/66_patients_riaz")
    merged_vcf = a.merge_vcdf()
    merged_vcf.to_csv("/Users/dch/Downloads/riaz_mutect2_66patients.tsv",sep='\t')
    #a = VCF_UTILS.parse_vcf_file("/Users/dch/Downloads/mutect2_riaz/annotated_Pt90.vcf","PT90")
    #print(a) 
