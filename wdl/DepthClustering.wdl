version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow ClusterDepth {
  input {
    File del_bed
    File dup_bed
    String batch
    File ploidy_table

    File contig_list
    File exclude_intervals
    Float exclude_overlap_fraction

    String? clustering_algorithm
    Float depth_interval_overlap
    Int depth_breakend_window

    File? gatk_to_svtk_script
    File? cnv_bed_to_gatk_vcf_script
    Float? java_mem_fraction

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_cnv_bed_to_gatk_vcf
    RuntimeAttr? runtime_override_concat_del_dup
    RuntimeAttr? runtime_attr_exclude_intervals_depth
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_attr_multi_gatk_to_svtk_vcf
  }

  call tasks_cluster.CNVBedToGatkVcf as DelBedToVcf {
    input:
      bed=del_bed,
      script=cnv_bed_to_gatk_vcf_script,
      contig_list=contig_list,
      ploidy_table=ploidy_table,
      output_prefix="~{batch}.depth.gatk_formatted.del",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_cnv_bed_to_gatk_vcf
  }

  call tasks_cluster.CNVBedToGatkVcf as DupBedToVcf {
    input:
      bed=dup_bed,
      script=cnv_bed_to_gatk_vcf_script,
      contig_list=contig_list,
      ploidy_table=ploidy_table,
      output_prefix="~{batch}.depth.gatk_formatted.dup",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_cnv_bed_to_gatk_vcf
  }

  call tasks_cohort.ConcatVcfs as ConcatDelDup {
    input:
      vcfs=[DelBedToVcf.out, DupBedToVcf.out],
      vcfs_idx=[DelBedToVcf.out_index, DupBedToVcf.out_index],
      allow_overlaps=true,
      outfile_prefix="~{batch}.depth.gatk_formatted.del_dup",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_del_dup
  }

  call tasks_cluster.ExcludeIntervalsDepth {
    input:
      vcf=ConcatDelDup.concat_vcf,
      overlap_fraction=exclude_overlap_fraction,
      output_prefix="~{batch}.depth.exclude_intervals",
      intervals=exclude_intervals,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_exclude_intervals_depth
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=[ExcludeIntervalsDepth.out],
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.depth.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=clustering_algorithm,
        depth_sample_overlap=0,
        depth_interval_overlap=depth_interval_overlap,
        depth_breakend_window=depth_breakend_window,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{batch}_depth_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
  }

  call tasks_cohort.ConcatVcfs as ConcatContigs {
    input:
      vcfs=SVCluster.out,
      vcfs_idx=SVCluster.out_index,
      naive=true,
      outfile_prefix="~{batch}.depth.clustered.gatk_formatted",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_del_dup
  }

  call tasks_cluster.GatkToSvtkVcf {
    input:
      vcf=ConcatContigs.concat_vcf,
      output_prefix="~{batch}.depth.clustered.svtk_formatted",
      script=gatk_to_svtk_script,
      source="depth",
      contig_list=contig_list,
      remove_formats="CN",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_multi_gatk_to_svtk_vcf
  }

  output {
    File clustered_vcf = ConcatContigs.concat_vcf
    File clustered_vcf_index = ConcatContigs.concat_vcf_idx
  }
}

