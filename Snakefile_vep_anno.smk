import re
import os
import glob

fasta_file="Data/project_data/raw/genome_files/hs37d5_spikein.fa"
vep_config="vep_GRCh37.config"
vep_cache_dir="<vep_cache_path>/vep/99/cachedir"
vep_conda_env="<conda_env_path>/vep/99"

def getAllSamples():
    all = glob.glob("Data/raw/vcfs*/vcfs*/*gz")
    ans = [re.sub(".vcf.gz", "", os.path.basename(x)) for x in all]
    return(ans)

def getInputVCF(wildcards):
    import os.path
    folders=["complex_care", "genomeClinic", "rna_vc"]
    for x in folders:
        out="Data/raw/vcfs_" + x + "/vcfs_" + x + "/" + wildcards['sampleID'] + ".vcf.gz"
        if os.path.isfile(out):
            return(out)

    raise Error("can not find VCF for sample: " + wildcards['sampleID'])
    

allSampleIDs = getAllSamples()


rule all:
    input: "Data/processed/processed_results/combined_rna_results_cagi6_with_variants_browseable.tsv.gz"
    output: touch("anno.done")


rule vep_anno:
    params:
        VEP_CONDA_ENV=vep_conda_env,
        VEP_CACHE_DIR=vep_cache_dir,
        BUFFER_SIZE=200000
    input: 
        VCF=getInputVCF,
        CONFIG=vep_config,
        FASTA=fasta_file
    output: 
        VEP="Data/processed/vep_annotation/{sampleID}.vep.vcf.gz"
    threads: 40
    resources:
        mem_mb=lambda wildcards, attempt, threads: threads * attempt * 800
    shell: """
        export PERL5LIB="{PARAM.VEP_CACHE_DIR}/Plugins/99_GRCh37/loftee:$PERL5LIB"
        export PERL5LIB="{PARAM.VEP_CAHCE_DIR}/Plugins/99_GRCh37/UTRannotator:$PERL5LIB"

        . <conda_install_path>/etc/profile.d/conda.sh
        
        conda activate {params.VEP_CONDA_ENV}
        src/bash/correct_vcf.sh {input.VCF} | \
            bcftools norm -m -both -f {input.FASTA} --check-ref ws --force --threads {threads} | \
            vep --config {input.CONFIG} \
                --dir {params.VEP_CACHE_DIR} \
                --fork {threads} \
                --buffer_size {params.BUFFER_SIZE} \
                --format vcf \
                --vcf \
                --minimal \
                --output_file STDOUT \
                --plugin SpliceAI,snv={Splice_AI_raw_data_folder}/spliceai_scores.raw.snv.hg19.vcf.gz,indel={Splice_AI_raw_data_folder}/spliceai_scores.raw.indel.hg19.vcf.gz | \
            bgzip -c -l 9 -@ {threads} > {output.VEP}
    """

rule normalize_vcf:
    input:
        VCF="{file}.vep.vcf.gz",
        TBI="{file}.vep.vcf.gz.tbi"
    output:
        NORM_VCF="{file}.norm.vcf.gz",
        NORM_TBI="{file}.norm.vcf.gz.tbi"
    params:
        REF = fasta_file
    shell:
        """
        bcftools norm -m-both {input.VCF} | grep -v -w "*" | \
        bcftools norm -f {params.REF} | bgzip -c > {output.NORM_VCF}
        tabix -f -p vcf {output.NORM_VCF} 
        """

rule extract_tsv_from_vep:
    input:  
        VCF="{file}.vep.norm.vcf.gz",
        TBI="{file}.vep.norm.vcf.gz.tbi",
        script="src/r/parse_vcf_files.R"
    output: TSV="{file}.rare_vars.tsv.gz"
    conda: "envs/r_variantanno.yaml"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt, threads: threads * attempt * 3000
    script: "src/r/parse_vcf_files.R"
    

rule merge_exp_with_variants:
    input: 
        EXP="Data/processed/processed_results/combined_rna_results_cagi6.tsv",
        VAR=expand("Data/processed/vep_annotation/{sampleID}.rare_vars.tsv.gz", sampleID=allSampleIDs),
        EVE="Data/processed/EVE_dt_all.Rds",
        Script="src/r/merge_variants_with_outliers.R"
    output:
        ALL   ="Data/processed/processed_results/combined_rna_results_cagi6_with_variants_all.tsv.gz",
        BROWSE="Data/processed/processed_results/combined_rna_results_cagi6_with_variants_browseable.tsv.gz"
    conda: "envs/r_parallel_tables.yaml"
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt, threads: threads * attempt * 5000
    script: "src/r/merge_variants_with_outliers.R"


