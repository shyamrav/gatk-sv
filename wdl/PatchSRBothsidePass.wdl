version 1.0

import "Structs.wdl"

workflow PatchSRBothsidePass {
    input {
        Array[File] batch_vcfs
        File cohort_vcf
        File updated_bothside_pass_list
        String cohort_name

        File patch_script

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_get_non_ref_vids
        RuntimeAttr? runtime_attr_calculate_support_frac
    }

    scatter (i in range(length(batch_vcfs))) {
        call GetNonRefVariantLists {
            input:
                batch_vcf=batch_vcfs[i],
                cohort_vcf=cohort_vcf,
                prefix="~{cohort_name}.non_ref_variants.shard_~{i}",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_non_ref_vids
        }
    }

    Int num_batches = length(batch_vcfs)
    call RecalculateBothsideSupportFractions {
        input:
            patch_script=patch_script,
            non_ref_vid_lists=GetNonRefVariantLists.out,
            updated_bothside_pass_list=updated_bothside_pass_list,
            num_batches=num_batches,
            prefix="~{cohort_name}.sr_bothside_support.patched",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_calculate_support_frac
    }

    output {
        File out = RecalculateBothsideSupportFractions.out
    }
}

task GetNonRefVariantLists {
    input {
        File batch_vcf
        File cohort_vcf
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([batch_vcf, cohort_vcf], "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bcftools query -l ~{batch_vcf} > samples.list
        bcftools view --samples-file samples.list ~{cohort_vcf} \
            | bcftools view -G -i 'SUM(AC)>0||SUM(FORMAT/SR_GT)>0' \
            | bcftools query -f '%ID\n' \
            > ~{prefix}.list
    >>>
    output {
        File out = "~{prefix}.list"
    }
}


task RecalculateBothsideSupportFractions {
    input {
        File patch_script
        Array[File] non_ref_vid_lists
        File updated_bothside_pass_list
        Int num_batches
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(non_ref_vid_lists, "GB") + size(updated_bothside_pass_list, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        python ~{patch_script} \
            ~{write_lines(non_ref_vid_lists)} \
            ~{updated_bothside_pass_list} \
            ~{num_batches} \
            > ~{prefix}.txt
    >>>
    output {
        File out = "~{prefix}.txt"
    }
}