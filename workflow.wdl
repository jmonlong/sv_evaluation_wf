version 1.0

workflow sv_evaluation {
    meta {
    author: "Jean Monlong"
        email: "jean.monlong@gmail.com"
        description: "Evaluate SV calls against a truth set using truvari and/or sveval"
    }

    parameter_meta {
        CALLS_VCF: "VCF with SV calls"
        TRUTH_VCF: "VCF with SV truthset"
        REFERENCE_FASTA: "Genome reference FASTA file"
        CONF_REGIONS_BED: "Optional. BED file defining confident regions to evaluate on"
        SIMPREP_BED: "Optional. Simple repeat annotation to help match SVs in repeats"
        SAMPLE: "Optional. Sample name. Should match between calls and truthset VCF. If empty, it will use the first sample in CALLS_VCF."
    }

    input {
        File CALLS_VCF
        File TRUTH_VCF
        File REFERENCE_FASTA
        File? CONF_REGIONS_BED
        File? SIMPREP_BED
        String SAMPLE = ""
    }

    call indexVcfFasta {
        input:
        call_vcf=CALLS_VCF,
        truth_vcf=TRUTH_VCF,
        reference_fa=REFERENCE_FASTA
    }
    
    call truvariEvaluate {
        input:
        call_vcf=CALLS_VCF,
        call_vcf_index=indexVcfFasta.call_vcf_index,
        truth_vcf=TRUTH_VCF,
        truth_vcf_index=indexVcfFasta.truth_vcf_index,
        reference_fa=REFERENCE_FASTA,
        reference_fa_index=indexVcfFasta.reference_index,
        confident_bed=CONF_REGIONS_BED
    }

    call svevalEvaluate {
        input:
        call_vcf=CALLS_VCF,
        call_vcf_index=indexVcfFasta.call_vcf_index,
        truth_vcf=TRUTH_VCF,
        truth_vcf_index=indexVcfFasta.truth_vcf_index,
        confident_bed=CONF_REGIONS_BED,
        reference_fa=REFERENCE_FASTA,
        reference_fa_index=indexVcfFasta.reference_index,
        simprep_bed=SIMPREP_BED,
        sample=SAMPLE
    }

    output {
        File truvariSummary = truvariEvaluate.summary
        File truvariTarball = truvariEvaluate.output_tarball
        File svevalSummary = svevalEvaluate.summary
        File svevalTarball = svevalEvaluate.output_tarball
    }
}



##
#### TASKS
##

task indexVcfFasta {
    input {
        File call_vcf
        File? truth_vcf
        File reference_fa
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="quay.io/biocontainers/samtools:1.16.1--h6899075_1"
        Int disk_size = 5 * round(size(call_vcf, 'G') + size(reference_fa, 'G')) + 20
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ln -s ~{call_vcf} call.vcf.gz
        ln -s ~{truth_vcf} truth.vcf.gz
        ln -s ~{reference_fa} ref.fa
        
        samtools faidx ref.fa
        tabix -p vcf call.vcf.gz
        tabix -p vcf truth.vcf.gz
    >>>
    output {
        File reference_index = "ref.fa.fai"
        File call_vcf_index = "call.vcf.gz.tbi"
        File truth_vcf_index = "truth.vcf.gz.tbi"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task truvariEvaluate {
    input {
        File call_vcf
        File call_vcf_index
        File truth_vcf
        File truth_vcf_index
        File reference_fa
        File reference_fa_index
        File? confident_bed
        String other_args = "-r 2000 -C 2000"
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="quay.io/biocontainers/truvari:3.5.0--pyhdfd78af_0"
        Int disk_size = 5 * round(size(call_vcf, 'G') + size(reference_fa, 'G')) + 20
    }
    Boolean confident_provided = defined(confident_bed)
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ln -s ~{call_vcf} call.vcf.gz
        ln -s ~{truth_vcf} truth.vcf.gz
        ln -s ~{call_vcf_index} call.vcf.gz.tbi
        ln -s ~{truth_vcf_index} truth.vcf.gz.tbi
        ln -s ~{reference_fa} ref.fa
        ln -s ~{reference_fa_index} ref.fa.fai

        ## confident regions if provided
        REG_ARGS=""
        if [ ~{confident_provided} == true ]
        then
        REG_ARGS="--includebed ~{confident_bed}"
        fi
           
        truvari bench --no-ref a -b truth.vcf.gz -c call.vcf.gz -f ref.fa ~{other_args} $REG_ARGS -o truvari_output/
        cp truvari_output/summary.txt summary.json
        tar -czvf truvari_output.tar.gz truvari_output
    >>>
    output {
        File summary = "summary.json"
        File output_tarball = "truvari_output.tar.gz"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task svevalEvaluate {
    input {
        File call_vcf
        File call_vcf_index
        File truth_vcf
        File truth_vcf_index
        File reference_fa
        File reference_fa_index
        File? confident_bed
        File? simprep_bed
        String sample = ""
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="quay.io/jmonlong/sveval:v2.3.0"
        Int disk_size = 5 * round(size(call_vcf, 'G') + size(reference_fa, 'G')) + 20
    }
    Boolean simprep_provided = defined(simprep_bed)
    Boolean confident_provided = defined(confident_bed)
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir vcf

        ## guess sample name if necessary
        SAMPLE=~{sample}
        if [ "$SAMPLE" == "" ]
        then
            SAMPLE=`bcftools view -h ~{call_vcf} | tail -1 | cut -f 10`
        fi

        ln -s ~{call_vcf} vcf/sveval-caller-${SAMPLE}.vcf.gz
        ln -s ~{call_vcf_index} vcf/sveval-caller-${SAMPLE}.vcf.gz.tbi
        
        ln -s ~{truth_vcf} vcf/sveval-truth-baseline.vcf.gz
        ln -s ~{truth_vcf_index} vcf/sveval-truth-baseline.vcf.gz.tbi
        
        ln -s ~{reference_fa} ref.fa
        ln -s ~{reference_fa_index} ref.fa.fai
        
        rm -f config.yaml
        echo 'ref_fa: "ref.fa"' >> config.yaml
        echo 'exp: "sveval"' >> config.yaml
        echo 'methods: "caller"' >> config.yaml
        echo 'out_prefix: "sveval_results"' >> config.yaml
        echo 'eval: "call geno"' >> config.yaml
        echo "samples: \"${SAMPLE}\"" >> config.yaml

        ## confident regions if provided
        if [ ~{confident_provided} == true ]
        then
            echo 'regions: "all conf"' >> config.yaml
            echo 'conf: "~{confident_bed}"' >> config.yaml
        else
            echo 'regions: "all"' >> config.yaml            
        fi

        ## simple repeat track if provided
        if [ ~{simprep_provided} == true ]
        then
            echo 'simprep_bed: "~{simprep_bed}"' >> config.yaml
        fi

        # ## REMOVE WHEN SNAKEMAKE HAS DEFAULTS
        # echo 'check_inv: False' >> config.yaml
        # echo 'min_cov: 0.5' >> config.yaml
        # echo -e "envm:\n  bgzip: 'bgzip'\n  bcftools: 'bcftools'\n  sveval: 'sveval'" >> config.yaml
        
        cp /sveval/snakemake/Snakefile .
        snakemake --configfile config.yaml --cores ~{threadCount}
        
        tar -czvf sveval_output.tar.gz sveval_results.pdf sveval_results-prcurve.tsv sveval_results-persize.tsv out/sveval*.RData
    >>>
    output {
        File summary = "sveval_results-prcurve.tsv"
        File output_tarball = "sveval_output.tar.gz"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}
