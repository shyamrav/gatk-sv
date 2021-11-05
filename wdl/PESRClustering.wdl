version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow ClusterPESR {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    File ploidy_table
    String batch
    String caller

    File? svtk_to_gatk_script

    File exclude_intervals
    File contigs

    Float pesr_interval_overlap
    Float pesr_breakend_window
    String? algorithm

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    Float? java_mem_fraction

    String gatk_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_override_concat_vcfs_pesr
  }

  call tasks_cluster.SvtkToGatkVcf {
    input:
      vcfs=vcfs,
      ploidy_table=ploidy_table,
      output_prefix="~{batch}.~{caller}.~{contig}.reformatted",
      script=svtk_to_gatk_script,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_exclude_intervals_pesr
  }

  call tasks_cluster.ExcludeIntervalsPESR {
    input:
      vcfs=SvtkToGatkVcf.out,
      output_prefix="~{batch}.~{caller}.exclude_intervals",
      intervals=exclude_intervals,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_exclude_intervals_pesr
  }

  Array[Array[String]] contiglist = read_tsv(contigs)
  scatter (contig in contiglist) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=ExcludeIntervalsPESR.out,
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.~{caller}.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=algorithm,
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

  output {
    File clustered_vcf = ConcatVcfs.concat_vcf
    File clustered_vcf_index = SVCluster.concat_vcf_idx
  }
}