version 1.0

task SVCluster {
    input {
        Array[File] vcfs
        File ploidy_table
        String output_prefix

        String? contig

        Boolean? fast_mode
        Boolean? omit_members
        Boolean? enable_cnv
        Boolean? default_no_call

        String? algorithm
        String? insertion_length_summary_strategy
        String? breakpoint_summary_strategy
        String? alt_allele_summary_strategy

        Float? defrag_padding_fraction
        Float? defrag_sample_overlap

        Float? depth_sample_overlap
        Float? depth_interval_overlap
        Int? depth_breakend_window
        Float? mixed_sample_overlap
        Float? mixed_interval_overlap
        Int? mixed_breakend_window
        Float? pesr_sample_overlap
        Float? pesr_interval_overlap
        Int? pesr_breakend_window

        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        Float? java_mem_fraction
        String? variant_prefix

        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcfs: {
                  localization_optional: true
              }
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcfs, "GB") * 2 + size(reference_fasta, "GB")),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail

        function getJavaMem() {
            # get JVM memory in GiB by getting total memory from /proc/meminfo
            # and multiplying by java_mem_fraction
            cat /proc/meminfo \
                | awk -v MEM_FIELD="$1" '{
                    f[substr($1, 1, length($1)-1)] = $2
                } END {
                    printf "%.2fG", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1048576
                }'
        }
        JVM_MAX_MEM=$(getJavaMem MemTotal)
        echo "JVM memory: $JVM_MAX_MEM"

        awk '{print "-V "$0}' ~{write_lines(vcfs)} > arguments.txt

        gatk --java-options "-Xmx${JVM_MAX_MEM}" SVCluster \
            --arguments_file arguments.txt \
            --output ~{output_prefix}.vcf.gz \
            --ploidy-table ~{ploidy_table} \
            --variant-prefix ~{variant_prefix} \
            --reference ~{reference_fasta} \
            ~{"-L " + contig} \
            ~{true="--fast-mode" false="" fast_mode} \
            ~{true="--enable-cnv" false="" enable_cnv} \
            ~{true="--omit-members" false="" omit_members} \
            ~{true="--default-no-call" false="" default_no_call} \
            ~{"--algorithm " + algorithm} \
            ~{"--defrag-padding-fraction " + defrag_padding_fraction} \
            ~{"--defrag-sample-overlap " + defrag_sample_overlap} \
            ~{"--depth-sample-overlap " + depth_sample_overlap} \
            ~{"--depth-interval-overlap " + depth_interval_overlap} \
            ~{"--depth-breakend-window " + depth_breakend_window} \
            ~{"--mixed-sample-overlap " + mixed_sample_overlap} \
            ~{"--mixed-interval-overlap " + mixed_interval_overlap} \
            ~{"--mixed-breakend-window " + mixed_breakend_window} \
            ~{"--pesr-sample-overlap " + pesr_sample_overlap} \
            ~{"--pesr-interval-overlap " + pesr_interval_overlap} \
            ~{"--pesr-breakend-window " + pesr_breakend_window} \
            ~{"--insertion-length-summary-strategy " + insertion_length_summary_strategy} \
            ~{"--breakpoint-summary-strategy " + breakpoint_summary_strategy} \
            ~{"--alt-allele-summary-strategy " + alt_allele_summary_strategy}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExcludeIntervalsPESR {
    input {
        Array[File] vcfs
        File intervals
        File? intervals_index
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcfs, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        Array[File] out = glob("out/~{output_prefix}.*.vcf.gz")
    }
    command <<<
        set -euo pipefail
        mkdir out/
        i=0
        while read vcf; do
            bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' $vcf \
                | awk '$1!="."' \
                > ends.bed
            bedtools intersect -wa -a ends.bed -b ~{intervals} | cut -f4 | sort | uniq | cut -f2 \
                > excluded_vids.list
            bcftools view -i '%ID!=@excluded_vids.list' $VCF -Oz -o out/~{output_prefix}.${i}.vcf.gz
            i=$((i+1))
        done < ~{write_lines(vcfs)}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SvtkToGatkVcf {
    input {
        Array[File] vcfs
        File ploidy_table
        File? script
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcfs, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        Array[File] out = glob("out/~{output_prefix}.*.vcf.gz")
    }
    command <<<
        set -euo pipefail
        mkdir out/
        i=0
        while read vcf; do
            python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
                $vcf \
                out/~{output_prefix}.$i.vcf.gz \
                ~{ploidy_table}
            i=$((i+1))
        done < ~{write_lines(vcfs)}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



task GatkToSvtkVcf {
    input {
        File vcf
        File? script
        String algorithm
        File contig_list
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
    }
    command <<<
        set -euo pipefail
        python ~{default="/opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py" script} \
            ~{vcf} \
            ~{output_prefix}.vcf.gz \
            ~{algorithm} \
            ~{contig_list}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}