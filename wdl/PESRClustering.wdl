version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow ClusterPESR {
  input {
    Array[File] vcfs

    File ploidy_table
    String batch
    String caller

    File? svtk_to_gatk_script
    File? gatk_to_svtk_script
    Float? java_mem_fraction

    File exclude_intervals
    File contig_list

    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? clustering_algorithm

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_multi_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_multi_gatk_to_svtk_vcf
  }

  call tasks_cluster.MultiSvtkToGatkVcf {
    input:
      vcfs=vcfs,
      ploidy_table=ploidy_table,
      output_suffix="gatk_formatted",
      script=svtk_to_gatk_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_multi_svtk_to_gatk_vcf
  }

  call tasks_cluster.MultiExcludeIntervalsPESR {
    input:
      vcfs=MultiSvtkToGatkVcf.out,
      output_suffix="exclude_intervals",
      intervals=exclude_intervals,
      intervals_index=exclude_intervals + ".tbi",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_exclude_intervals_pesr
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=MultiExcludeIntervalsPESR.out,
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.~{caller}.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=clustering_algorithm,
        pesr_sample_overlap=0,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{batch}_~{caller}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVCluster.out,
      vcfs_idx=SVCluster.out_index,
      naive=true,
      outfile_prefix="~{batch}.~{caller}.clustered",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs_pesr
  }

  call tasks_cluster.GatkToSvtkVcf {
    input:
      vcf=ConcatVcfs.concat_vcf,
      output_prefix="~{batch}.~{caller}.clustered.svtk_formatted",
      script=gatk_to_svtk_script,
      source=caller,
      contig_list=contig_list,
      remove_formats="CN",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_multi_gatk_to_svtk_vcf
  }

  output {
    File clustered_vcf = GatkToSvtkVcf.out
    File clustered_vcf_index = GatkToSvtkVcf.out_index
  }
}