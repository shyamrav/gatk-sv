## This WDL pipeline implements SV calling of mobile elements with Scramble by Rebecca I. Torene.
## https://github.com/GeneDx/scramble
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format

version 1.0

import "Structs.wdl"

workflow Scramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File reference_fasta
    Boolean detect_deletions
    String scramble_docker
    RuntimeAttr? runtime_attr_scramble
  }
    
  parameter_meta {
    bam_or_cram_file: "A .bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "Index for bam_or_cram_file."
    sample_name: "A sample name. Outputs will be sample_name+'.scramble.insertions.vcf.gz' and sample_name+'.scramble.deletions.vcf.gz'."
    reference_fasta: "A FASTA file with the reference used to align bam or cram file."
    detect_deletions: "Run deletion detection as well as mobile element insertion."
  }
  
  meta {
    author: "Ted Sharpe, et al"
    email: "tsharpe@broadinstitute.org"
  }

  call RunScramble {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      detect_deletions = detect_deletions,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble
  }

  output {
    File insertions_vcf = RunScramble.insertions_vcf
    File insertions_index = RunScramble.insertions_index
    File? deletions_vcf = RunScramble.deletions_vcf
    File? deletions_index = RunScramble.deletions_index
  }
}

task RunScramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File reference_fasta
    Boolean detect_deletions
    String scramble_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File insertions_vcf = "~{sample_name}.scramble.insertions.vcf.gz"
    File insertions_index = "~{sample_name}.scramble.insertions.vcf.gz.tbi"
    File? deletions_vcf = "~{sample_name}.scramble.deletions.vcf.gz"
    File? deletions_index = "~{sample_name}.scramble.deletions.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    xDir=$PWD
    scrambleDir="/app"

    cat ~{reference_fasta} | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref
    $scrambleDir/cluster_identifier/src/build/cluster_identifier ~{bam_or_cram_file} > clusters.txt
    split -a3 -l3000 clusters.txt xyzzy

    for fil in xyzzy???
      do Rscript --vanilla $scrambleDir/cluster_analysis/bin/SCRAMble.R --out-name $xDir/$fil \
           --cluster-file $xDir/$fil --install-dir $scrambleDir/cluster_analysis/bin \
           --mei-refs $scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa \
           --ref $xDir/ref --no-vcf --eval-meis ~{true='--eval-dels' false='' detect_deletions}
    done

    echo '##fileformat=VCFv4.3' > hdr
    echo '##reference='"~{reference_fasta}" >> hdr
    echo '##source=scramble' >> hdr
    echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">' >> hdr
    echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> hdr
    echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' >> hdr
    echo '##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">' >> hdr
    echo '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness [++,+-,-+,--]">' >> hdr
    echo '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">' >> hdr
    echo '##INFO=<ID=MEI_START,Number=1,Type=Integer,Description="Start of alignment to canonical MEI sequence">' >> hdr
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> hdr
    blastdbcmd -db ref -entry all -outfmt '##contig=<ID=%a,length=%l>' >> hdr
    echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'"~{sample_name}" >> hdr

    cp hdr insertsVCF.tmp

    echo '{ FS=OFS="\t" }' > processInserts.awk
    echo '{ if(FNR<2)next' >> processInserts.awk
    echo '  split($1,loc,":")' >> processInserts.awk
    echo '  start=loc[2]+1' >> processInserts.awk
    echo '  len=length($8)' >> processInserts.awk
    echo '  end=start+len-1' >> processInserts.awk
    echo '  print loc[1],start,".","N","<INS:ME:" toupper($2) ">",int($6),"PASS","END=" end ";SVTYPE=INS;SVLEN=" len ";MEI_START=" $10 ";STRANDS=+-;CHR2=" loc[1] ";ALGORITHMS=scramble","GT","0/1" }' >> processInserts.awk

    awk -f processInserts.awk xyzzy???_MEIs.txt >> insertsVCF.tmp
    bcftools sort -Oz <insertsVCF.tmp >"~{sample_name}.scramble.insertions.vcf.gz"
    bcftools index -ft "~{sample_name}.scramble.insertions.vcf.gz"

    if [ ~{detect_deletions} == "true" ]
    then
      cp hdr deletesVCF.tmp

      echo '{ FS=OFS="\t" }' > processDeletes.awk
      echo '{ if(FNR<2)next' >> processDeletes.awk
      echo '  Q= $11=="NA" ? ($15=="NA"?".":$15) : ($15=="NA"?$11:($11+$15)/2)' >> processDeletes.awk
      echo '  print $1,$2+1,".","N","<DEL>",Q=="."?".":int(Q),"PASS","END=" $3+1 ";SVTYPE=DEL;SVLEN=" $5 ";STRANDS=+-;CHR2=" $1 ";ALGORITHMS=scramble","GT","0/1" }' >> processDeletes.awk

      awk -f processDeletes.awk xyzzy???_PredictedDeletions.txt >> deletesVCF.tmp
      bcftools sort -Oz <deletesVCF.tmp >"~{sample_name}.scramble.deletions.vcf.gz"
      bcftools index -ft "~{sample_name}.scramble.deletions.vcf.gz"
    fi
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: scramble_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
