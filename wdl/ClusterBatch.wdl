version 1.0

import "PESRClustering.wdl" as pesr
import "DepthClustering.wdl" as depth
import "ClusterBatchMetrics.wdl" as metrics

workflow ClusterBatch {
  input {
    Array[File]? manta_vcfs
    Array[File]? wham_vcfs
    Array[File]? melt_vcfs
    File del_bed
    File dup_bed
    String batch
    File ploidy_table

    File contig_list

    File? cnv_bed_to_gatk_vcf_script
    File? svtk_to_gatk_script
    File? gatk_to_svtk_script

    File depth_exclude_intervals
    Float depth_exclude_overlap_fraction
    String? depth_clustering_algorithm
    Float depth_interval_overlap
    Int depth_breakend_window

    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    Float? java_mem_fraction
    String batch

    Int pesr_svsize
    Float pesr_frac
    String pesr_flags
    Int pesr_distance
    File pesr_exclude_list

    File pesr_exclude_list
    File? depth_exclude_list
    Float? depth_exclude_list_frac_max

    String depth_flags
    Float depth_frac

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? linux_docker  # required if run_module_metrics = true
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true
    File? baseline_depth_vcf  # baseline files are optional for metrics workflow
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_multi_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_attr_svcluster_pesr
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_multi_gatk_to_svtk_vcf
    RuntimeAttr? runtime_attr_cnv_bed_to_gatk_vcf
    RuntimeAttr? runtime_override_concat_del_dup
    RuntimeAttr? runtime_attr_exclude_intervals_depth
    RuntimeAttr? runtime_attr_svcluster_depth
  }

  if (defined(manta_vcfs) && (length(select_first([manta_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_manta {
      input:
        vcfs=select_first([manta_vcfs]),
        ploidy_table=ploidy_table,
        batch=batch,
        caller="manta",
        svtk_to_gatk_script=svtk_to_gatk_script,
        gatk_to_svtk_script=gatk_to_svtk_script,
        exclude_intervals=pesr_exclude_list,
        contig_list=contig_list,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        clustering_algorithm=pesr_clustering_algorithm,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_multi_svtk_to_gatk_vcf=runtime_attr_multi_svtk_to_gatk_vcf,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr,
        runtime_attr_svcluster=runtime_attr_svcluster_pesr,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_multi_gatk_to_svtk_vcf=runtime_attr_multi_gatk_to_svtk_vcf
    }
  }

  if (defined(wham_vcfs) && (length(select_first([wham_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_wham {
      input:
        vcfs=select_first([wham_vcfs]),
        ploidy_table=ploidy_table,
        batch=batch,
        caller="wham",
        svtk_to_gatk_script=svtk_to_gatk_script,
        gatk_to_svtk_script=gatk_to_svtk_script,
        exclude_intervals=pesr_exclude_list,
        contig_list=contig_list,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        clustering_algorithm=pesr_clustering_algorithm,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_multi_svtk_to_gatk_vcf=runtime_attr_multi_svtk_to_gatk_vcf,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr,
        runtime_attr_svcluster=runtime_attr_svcluster_pesr,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_multi_gatk_to_svtk_vcf=runtime_attr_multi_gatk_to_svtk_vcf
    }
  }

  if (defined(melt_vcfs) && (length(select_first([melt_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_melt {
      input:
        vcfs=select_first([melt_vcfs]),
        ploidy_table=ploidy_table,
        batch=batch,
        caller="melt",
        svtk_to_gatk_script=svtk_to_gatk_script,
        gatk_to_svtk_script=gatk_to_svtk_script,
        exclude_intervals=pesr_exclude_list,
        contig_list=contig_list,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        clustering_algorithm=pesr_clustering_algorithm,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_multi_svtk_to_gatk_vcf=runtime_attr_multi_svtk_to_gatk_vcf,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr,
        runtime_attr_svcluster=runtime_attr_svcluster_pesr,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_multi_gatk_to_svtk_vcf=runtime_attr_multi_gatk_to_svtk_vcf
    }
  }

  call depth.ClusterDepth as ClusterDepth {
  	input:
      del_bed=del_bed,
      dup_bed=dup_bed,
      batch=batch,
      ploidy_table=ploidy_table,
      contig_list=contig_list,
      exclude_intervals=depth_exclude_intervals,
      exclude_overlap_fraction=depth_exclude_overlap_fraction,
      clustering_algorithm=depth_clustering_algorithm,
      depth_interval_overlap=depth_interval_overlap,
      depth_breakend_window=depth_breakend_window,
      gatk_to_svtk_script=gatk_to_svtk_script,
      cnv_bed_to_gatk_vcf_script=cnv_bed_to_gatk_vcf_script,
      java_mem_fraction=java_mem_fraction,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_cnv_bed_to_gatk_vcf=runtime_attr_cnv_bed_to_gatk_vcf,
      runtime_override_concat_del_dup=runtime_override_concat_del_dup,
      runtime_attr_exclude_intervals_depth=runtime_attr_exclude_intervals_depth,
      runtime_attr_svcluster=runtime_attr_svcluster_depth,
      runtime_attr_multi_gatk_to_svtk_vcf=runtime_attr_multi_gatk_to_svtk_vcf
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call metrics.ClusterBatchMetrics {
      input:
        name = batch,
        depth_vcf = ClusterDepth.clustered_vcf,
        manta_vcf = ClusterPESR_manta.clustered_vcf,
        wham_vcf = ClusterPESR_wham.clustered_vcf,
        melt_vcf = ClusterPESR_melt.clustered_vcf,
        baseline_depth_vcf = baseline_depth_vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        contig_list = select_first([primary_contigs_list]),
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker]),
        linux_docker = select_first([linux_docker])
    }
  }

  output {
    File depth_vcf = ClusterDepth.clustered_vcf
    File? manta_vcf = ClusterPESR_manta.clustered_vcf
    File? wham_vcf = ClusterPESR_wham.clustered_vcf
    File? melt_vcf = ClusterPESR_melt.clustered_vcf

    File? metrics_file_clusterbatch = ClusterBatchMetrics.metrics_file
  }
}
