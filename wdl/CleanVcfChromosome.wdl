version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "CleanVcf5.wdl" as c5
import "HailMerge.wdl" as HailMerge

workflow CleanVcfChromosome {
  input {
    File vcf
    String contig
    File background_list
    File ped_file
    File allosome_fai
    String prefix
    Int max_shards_per_chrom_step1
    File bothsides_pass_list
    Int min_records_per_shard_step1
    Int samples_per_step2_shard
    Int clean_vcf5_records_per_shard
    Int? clean_vcf5_threads_per_task
    File? outlier_samples_list
    Int? max_samples_per_shard_step3

    File hail_script
    String project

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_updates_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_1b
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5_scatter
    RuntimeAttr? runtime_override_clean_vcf_5_make_cleangq
    RuntimeAttr? runtime_override_clean_vcf_5_find_redundant_multiallelics
    RuntimeAttr? runtime_override_clean_vcf_5_polish
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup

    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_fix_header

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4
    RuntimeAttr? runtime_override_clean_vcf_5_scatter
    RuntimeAttr? runtime_override_clean_vcf_5_make_cleangq
    RuntimeAttr? runtime_override_clean_vcf_5_find_redundant_multiallelics
    RuntimeAttr? runtime_override_clean_vcf_5_polish

  }

  call MiniTasks.SplitVcf as SplitVcfToClean {
    input:
      vcf=vcf,
      contig=contig,
      prefix="~{prefix}.~{contig}.shard_",
      n_shards=max_shards_per_chrom_step1,
      min_vars_per_shard=min_records_per_shard_step1,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_vcf_to_clean
  }

  scatter ( vcf_shard in SplitVcfToClean.vcf_shards ) {
    call CleanVcf1a {
      input:
        vcf=vcf_shard,
        background_fail_list=background_list,
        bothsides_pass_list=bothsides_pass_list,
        ped_file=ped_file,
        allosome_fai=allosome_fai,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_1a
    }
  }

  call HailMerge.HailMerge as CombineStep1Vcfs {
    input:
      vcfs=CleanVcf1a.intermediate_vcf,
      prefix="~{prefix}.cleanVCF_step1.intermediate_vcf.merged",
      hail_script=hail_script,
      project=project,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_override_preconcat=runtime_override_preconcat,
      runtime_override_hail_merge=runtime_override_hail_merge,
      runtime_override_fix_header=runtime_override_fix_header
  }

  call MiniTasks.CatUncompressedFiles as CombineStep1SexChrRevisions {
    input:
      shards=CleanVcf1a.sex,
      outfile_name=prefix + ".cleanVCF_step1.sexchr_revise.merged.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_step_1_sex_chr_revisions
  }

  call CleanVcf1b {
    input:
      intermediate_vcf=CombineStep1Vcfs.merged_vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_1b
  }

  call MiniTasks.SplitUncompressed as SplitIncludeList {
    input:
      whole_file=CleanVcf1a.include_list[0],
      lines_per_shard=samples_per_step2_shard,
      shard_prefix="includeexclude.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_include_list
  }

  scatter ( included_interval in SplitIncludeList.shards ){
    call CleanVcf2 {
      input:
        normal_revise_vcf=CleanVcf1b.normal,
        include_list=included_interval,
        multi_cnvs=CleanVcf1b.multi,
        vcftools_idx=CleanVcf1b.vcftools_idx,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_2
      }
  }

  call MiniTasks.CatUncompressedFiles as CombineCleanVcf2 {
    input:
      shards=CleanVcf2.out,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_clean_vcf_2
  }

  call CleanVcf3 {
    input:
      rd_cn_revise=CombineCleanVcf2.outfile,
      max_samples_shard = max_samples_per_shard_step3,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_3
  }

  scatter ( rd_cn_revise in CleanVcf3.shards ){
    call CleanVcf4 {
      input:
        rd_cn_revise=rd_cn_revise,
        normal_revise_vcf=CleanVcf1b.normal,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_4
    }
  }

  call MiniTasks.CatUncompressedFiles as CombineRevised4 {
    input:
      shards=CleanVcf4.out,
      outfile_name="revise.vcf.lines.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_revised_4
  }

  call MiniTasks.CatUncompressedFiles as CombineMultiIds4 {
    input:
      shards=CleanVcf4.multi_ids,
      outfile_name="multi.geno.ids.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_multi_ids_4
  }

  call c5.CleanVcf5 {
    input:
      revise_vcf_lines=CombineRevised4.outfile,
      normal_revise_vcf=CleanVcf1b.normal,
      ped_file=ped_file,
      sex_chr_revise=CombineStep1SexChrRevisions.outfile,
      multi_ids=CombineMultiIds4.outfile,
      outlier_samples_list=outlier_samples_list,
      contig=contig,
      prefix=prefix,
      records_per_shard=clean_vcf5_records_per_shard,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override_scatter=runtime_override_clean_vcf_5_scatter,
      runtime_attr_override_make_cleangq=runtime_override_clean_vcf_5_make_cleangq,
      runtime_attr_override_find_redundant_multiallelics=runtime_override_clean_vcf_5_find_redundant_multiallelics,
      runtime_attr_override_polish=runtime_override_clean_vcf_5_polish
  }

  call DropRedundantCnvs {
    input:
      vcf=CleanVcf5.polished,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_updates_docker
  }

  call StitchFragmentedCnvs {
    input:
      vcf=DropRedundantCnvs.cleaned_vcf_shard,
      prefix="~{prefix}.stitched",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_stitch_fragmented_cnvs
  }

  call FinalCleanup {
    input:
      vcf=StitchFragmentedCnvs.stitched_vcf_shard,
      contig=contig,
      prefix="~{prefix}.final",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_final_cleanup

  }
  
  output {
    File out=FinalCleanup.final_cleaned_shard
    File out_idx=FinalCleanup.final_cleaned_shard_idx
  }
}


task CleanVcf1a {
  input {
    File vcf
    File background_fail_list
    File bothsides_pass_list
    File ped_file
    File allosome_fai
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, background_fail_list, bothsides_pass_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 2),
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

    # outputs
    # includelist.txt: the names of all the samples in the input vcf
    # sexchr.revise.txt: the names of the events where genotypes got tweaked on allosomes
    # int.vcf.gz: a revised vcf, bgzipped
    zcat ~{vcf} \
      | awk -v allosomeFile="~{allosome_fai}" -v pedFile="~{ped_file}" -v bgdFile="~{background_fail_list}" \
        -f /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1.awk \
      | bgzip \
      > int.vcf.gz
    rm ~{vcf} ~{background_fail_list} ~{ped_file}

    /opt/sv-pipeline/04_variant_resolution/scripts/add_bothsides_support_filter.py \
      --bgzip \
      --outfile int.w_bothsides.vcf.gz \
      int.vcf.gz \
      ~{bothsides_pass_list}
    tabix int.w_bothsides.vcf.gz
  >>>

  output {
    File include_list="includelist.txt"
    File sex="sexchr.revise.txt"
    File intermediate_vcf="int.w_bothsides.vcf.gz"
    File intermediate_vcf_idx="int.w_bothsides.vcf.gz.tbi"
  }
}


task CleanVcf1b {
  input {
    File intermediate_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intermediate_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 3),
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
    python /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1b.py ~{intermediate_vcf}
  >>>

  output {
    File multi = "multi.cnvs.txt"
    File normal = "normal.revise.vcf.gz"
    File vcftools_idx = "normal.revise.vcf.gz.csi"
  }
}

task CleanVcf2 {
  input {
    File normal_revise_vcf
    File include_list
    File multi_cnvs
    File vcftools_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size([normal_revise_vcf, include_list, multi_cnvs, vcftools_idx], "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 4.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + input_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
    set -eu -o pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part2.sh \
      ~{normal_revise_vcf} \
      ~{include_list} \
      ~{multi_cnvs} \
      "output.txt"
  >>>

  output {
    File out="output.txt"
  }
}


task CleanVcf3 {
  input {
    File rd_cn_revise
    Int? max_samples_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  Int max_samples_shard_ = select_first([max_samples_shard, 6000])
  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(rd_cn_revise, "GB")
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
    
    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part3.py ~{rd_cn_revise} -s ~{max_samples_shard_}

    # Ensure there is at least one shard
    if [ -z "$(ls -A shards/)" ]; then
      touch shards/out.0_0.txt
    fi
  >>>

  output {
     Array[File] shards = glob("shards/*")
  }
}


task CleanVcf4 {
  input {
    File rd_cn_revise
    File normal_revise_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([rd_cn_revise, normal_revise_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0 + input_size * 3.0,
                                  disk_gb: 50,
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
    python3 <<CODE
    import pysam
    import os

    # Inputs
    REGENO_FILE="~{rd_cn_revise}"
    VCF_FILE="~{normal_revise_vcf}"

    # Build map of variants to regenotype
    with open(REGENO_FILE) as f:
      vid_sample_cn_map = {}
      for line in f:
        tokens = line.strip().split('\t')
        vid = tokens[0]
        if vid not in vid_sample_cn_map:
          vid_sample_cn_map[vid] = []
        vid_sample_cn_map[vid].append(tuple(tokens[1:]))

    # Traverse VCF and replace genotypes
    with open("revise.vcf.lines.txt", "w") as f:
      vcf = pysam.VariantFile(VCF_FILE)
      num_vcf_records = 0
      for record in vcf:
        num_vcf_records += 1
        if record.id not in vid_sample_cn_map:
          continue
        for entry in vid_sample_cn_map[record.id]:
          s = record.samples[entry[0]]
          s['GT'] = (0, 1)
          s['RD_CN'] = int(entry[1])
        f.write(str(record))
      vcf.close()

    # Get batch size
    regeno_file_name_tokens = os.path.basename(REGENO_FILE).split('.')[1].split('_')
    batch_num = max(int(regeno_file_name_tokens[0]), 1)
    total_batch = max(int(regeno_file_name_tokens[1]), 1)
    segments = num_vcf_records / float(total_batch)
    print("{} {} {}".format(batch_num, total_batch, segments))

    vcf = pysam.VariantFile(VCF_FILE)
    # Max sample count with PE or SR GT over 3
    max_vf = max(len(vcf.header.samples) * 0.01, 2)
    record_start = (batch_num - 1) * segments
    record_end = batch_num * segments
    record_idx = 0
    print("{} {} {}".format(max_vf, record_start, record_end))
    multi_geno_ids = set([])
    for record in vcf:
      record_idx += 1
      if record_idx < record_start:
        continue
      elif record_idx > record_end:
        break
      num_gt_over_2 = 0
      for sid in record.samples:
        s = record.samples[sid]
        # Pick best GT
        if s['PE_GT'] is None:
          continue
        elif s['SR_GT'] is None:
          gt = s['PE_GT']
        elif s['PE_GT'] > 0 and s['SR_GT'] == 0:
          gt = s['PE_GT']
        elif s['PE_GT'] == 0:
          gt = s['SR_GT']
        elif s['PE_GQ'] >= s['SR_GQ']:
          gt = s['PE_GT']
        else:
          gt = s['SR_GT']
        if gt > 2:
          num_gt_over_2 += 1
          if record.id == "gnomad-sv-v3-TEST-SMALL.chr22_BND_chr22_173":
            print("{} {}".format(sid, num_gt_over_2))
            print("{} {} {} {}".format(s['PE_GT'], s['PE_GQ'], s['SR_GT'], s['SR_GQ']))
      if num_gt_over_2 > max_vf:
        multi_geno_ids.add(record.id)
    vcf.close()

    multi_geno_ids = sorted(list(multi_geno_ids))
    with open("multi.geno.ids.txt", "w") as f:
      for vid in multi_geno_ids:
        f.write(vid + "\n")
    CODE

    bgzip revise.vcf.lines.txt
    gzip multi.geno.ids.txt
  >>>

  output {
    File out="revise.vcf.lines.txt.gz"
    File multi_ids="multi.geno.ids.txt.gz"
  }
}


# Remove CNVs that are redundant with CPX events or other CNVs
task DropRedundantCnvs {
  input {
    File vcf
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String concise_vcf_name = contig + ".shard.no_CNV_redundancies.vcf.gz"

  Float input_size = size(vcf, "GiB")
  # disk is cheap, read/write speed is proportional to disk size, and disk IO is a significant time factor:
  # in tests on large VCFs, memory usage is ~1.0 * input VCF size
  # the biggest disk usage is at the end of the task, with input + output VCF on disk
  Int cpu_cores = 2 # speed up compression / decompression of VCFs
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75 + input_size * 1.5,
    disk_gb: ceil(100.0 + input_size * 2.0),
    cpu_cores: cpu_cores,
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
    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
      ~{vcf} ~{concise_vcf_name} --temp-dir ./tmp
  >>>

  output {
    File cleaned_vcf_shard = concise_vcf_name
  }
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 2),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  Float mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  runtime {
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    echo "First pass..."
    java -Xmx~{java_mem_mb}M -jar /opt/sv-pipeline/java/StitchFragmentedCNVs.jar 0.2 200000 0.2 ~{vcf} \
      | bgzip \
      > tmp.vcf.gz
    rm ~{vcf}

    echo "Second pass..."
    java -Xmx~{java_mem_mb}M -jar /opt/sv-pipeline/java/StitchFragmentedCNVs.jar 0.2 200000 0.2 tmp.vcf.gz \
      | bgzip \
      > ~{prefix}.vcf.gz
  >>>

  output {
    File stitched_vcf_shard = "~{prefix}.vcf.gz"
  }
}


# Final VCF cleanup
task FinalCleanup {
  input {
    File vcf
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(vcf, "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + input_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{vcf} stdout \
      | fgrep -v "##INFO=<ID=HIGH_SR_BACKGROUND" \
      | /opt/sv-pipeline/04_variant_resolution/scripts/sanitize_filter_field.py stdin stdout \
      | fgrep -v "##INFO=<ID=MEMBERS,Number=.,Type=String," \
      | bgzip -c \
      > ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File final_cleaned_shard = "~{prefix}.vcf.gz"
    File final_cleaned_shard_idx = "~{prefix}.vcf.gz.tbi"
  }
}