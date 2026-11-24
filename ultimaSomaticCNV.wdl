version 1.0

# Adapted from the Ultima Genomics controlFREEC_pipeline workflow (healthomics-workflows v1.37.1,
# Apache License 2.0). Tools run from the upstream container images through Apptainer.
# Paths are rooted at environment variables set by the modules each task loads:
# HG38_ULTIMA_ROOT (hg38-ultima) for the reference, HG38_UGFREEC_RESOURCES_ROOT
# (hg38-ugfreec-resources) for the calling regions, windows, SNP sites, mappability and
# low-confidence regions, UGFREEC_ROOT (ugfreec) for the container images.

struct GenomeResources {
  String genomeModule
  String refFasta
  String refFai
  String chrLenFile
  String genomeWindows
  String intervalList
  String snpVcf
  String gemMappability
  String cnvLcrBed
}

workflow ultimaSomaticCNV {
  input {
    Array[File]+ tumorCrams
    Array[File]+ tumorCramIndexes
    Array[File]+ normalCrams
    Array[File]+ normalCramIndexes
    String outputFileNamePrefix
    String reference = "hg38"
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
    String? sex
    Boolean highSensitivityMode = false
    Int numShards = 100
    Int scatterIntervalsBreak = 10000000
    Int minMapq = 1
  }

  parameter_meta {
    tumorCrams: "Tumor CRAM or BAM files of one sample"
    tumorCramIndexes: "Indexes of tumorCrams, in the same order"
    normalCrams: "Matched-normal CRAM or BAM files of one sample"
    normalCramIndexes: "Indexes of normalCrams, in the same order"
    outputFileNamePrefix: "Prefix for output file names"
    reference: "Reference genome build, selects the genome resources: hg38"
    chromosomes: "Chromosomes to call, a subset of chr1-chr22 and chrX; restrict it to run a test on part of the genome"
    sex: "Sample sex, XX or XY, for the expected copy number of chrX; Control-FREEC infers it when unset"
    highSensitivityMode: "Call from segment fold changes with ploidy 4 and 10% contamination assumed, for low-purity tumors"
    numShards: "Number of genomic shards for the SNP pileups"
    scatterIntervalsBreak: "Break scatter intervals at multiples of this many bases"
    minMapq: "Minimum mapping quality of reads counted toward coverage and SNP pileups"
  }

  String freecImage = "$UGFREEC_ROOT/images/ultimagenomics_ugbio_freec_1.28.0.sif"
  String cnvImage = "$UGFREEC_ROOT/images/ultimagenomics_ugbio_cnv_1.33.0.sif"

  Map[String, GenomeResources] genomeResources = {
    "hg38": {
      "genomeModule": "hg38-ultima/v0 hg38-ugfreec-resources/1.37.1",
      "refFasta": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta",
      "refFai": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta.fai",
      "chrLenFile": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/Homo_sapiens_assembly38.chr1-23.genome",
      "genomeWindows": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/Homo_sapiens_assembly38.w1000.chr1-23.bed",
      "intervalList": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list",
      "snpVcf": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/af-only-gnomad.hg38.AF_gt0.35.CHR1-24.vcf.gz",
      "gemMappability": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/out100m2_hg38.gem",
      "cnvLcrBed": "$HG38_UGFREEC_RESOURCES_ROOT/hg38/UG-High-Confidence-Regions/v2.1/ug_cnv_lcr.bed"
    }
  }

  GenomeResources genome = genomeResources[reference]
  String tumorName = basename(tumorCrams[0])
  String normalName = basename(normalCrams[0])

  call prepareRegions {
    input:
      chromosomes = chromosomes,
      chrLenFile = genome.chrLenFile,
      genomeWindows = genome.genomeWindows,
      intervalList = genome.intervalList,
      genomeModule = genome.genomeModule
  }

  call scatterIntervals {
    input:
      intervalList = prepareRegions.selectedIntervalList,
      numShards = numShards,
      breakBandsAtMultiplesOf = scatterIntervalsBreak,
      image = freecImage
  }

  scatter (shardBed in scatterIntervals.shardBeds) {
    call createMpileup as tumorMpileup {
      input:
        crams = tumorCrams,
        cramIndexes = tumorCramIndexes,
        shardBed = shardBed,
        outputPrefix = tumorName + "." + basename(shardBed, ".bed"),
        minMapq = minMapq,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        snpVcf = genome.snpVcf,
        image = freecImage
    }

    call createMpileup as normalMpileup {
      input:
        crams = normalCrams,
        cramIndexes = normalCramIndexes,
        shardBed = shardBed,
        outputPrefix = normalName + "." + basename(shardBed, ".bed"),
        minMapq = minMapq,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        snpVcf = genome.snpVcf,
        image = freecImage
    }
  }

  call concatFiles as tumorConcatPileups {
    input:
      files = tumorMpileup.pileup,
      outputFileName = tumorName + "_minipileup.pileup"
  }

  call concatFiles as normalConcatPileups {
    input:
      files = normalMpileup.pileup,
      outputFileName = normalName + "_minipileup.pileup"
  }

  scatter (tumor in zip(tumorCrams, tumorCramIndexes)) {
    scatter (chromosome in prepareRegions.orderedChromosomes) {
      call collectCoverage as tumorCoverage {
        input:
          cram = tumor.left,
          cramIndex = tumor.right,
          chromosome = chromosome,
          minMapq = minMapq,
          genomeModule = genome.genomeModule,
          refFasta = genome.refFasta,
          image = freecImage
      }
    }

    call concatFiles as tumorConcatCoverage {
      input:
        files = tumorCoverage.bedGraph,
        outputFileName = basename(tumor.left) + ".bedGraph"
    }
  }

  scatter (normal in zip(normalCrams, normalCramIndexes)) {
    scatter (chromosome in prepareRegions.orderedChromosomes) {
      call collectCoverage as normalCoverage {
        input:
          cram = normal.left,
          cramIndex = normal.right,
          chromosome = chromosome,
          minMapq = minMapq,
          genomeModule = genome.genomeModule,
          refFasta = genome.refFasta,
          image = freecImage
      }
    }

    call concatFiles as normalConcatCoverage {
      input:
        files = normalCoverage.bedGraph,
        outputFileName = basename(normal.left) + ".bedGraph"
    }
  }

  call bedGraphToCpn as tumorCpn {
    input:
      bedGraphs = tumorConcatCoverage.mergedFile,
      genomeWindows = prepareRegions.selectedWindows,
      chrLenFile = prepareRegions.selectedChrLenFile,
      outputPrefix = tumorName,
      image = freecImage
  }

  call bedGraphToCpn as normalCpn {
    input:
      bedGraphs = normalConcatCoverage.mergedFile,
      genomeWindows = prepareRegions.selectedWindows,
      chrLenFile = prepareRegions.selectedChrLenFile,
      outputPrefix = normalName,
      image = freecImage
  }

  call runControlFreec {
    input:
      sampleName = outputFileNamePrefix,
      tumorCram = tumorCrams[0],
      tumorPileup = tumorConcatPileups.mergedFile,
      tumorCpn = tumorCpn.cpn,
      normalCram = normalCrams[0],
      normalPileup = normalConcatPileups.mergedFile,
      normalCpn = normalCpn.cpn,
      chrLenFile = prepareRegions.selectedChrLenFile,
      sex = sex,
      highSensitivityMode = highSensitivityMode,
      genomeModule = genome.genomeModule,
      refFasta = genome.refFasta,
      refFai = genome.refFai,
      snpVcf = genome.snpVcf,
      gemMappability = genome.gemMappability,
      image = freecImage
  }

  call annotateCnvs {
    input:
      sampleName = outputFileNamePrefix,
      cnvCalls = runControlFreec.tumorCnvs,
      segments = runControlFreec.segments,
      highSensitivityMode = highSensitivityMode,
      genomeModule = genome.genomeModule,
      refFai = genome.refFai,
      cnvLcrBed = genome.cnvLcrBed,
      image = cnvImage
  }

  call filterCnvs {
    input:
      sampleName = outputFileNamePrefix,
      annotatedCnvBed = annotateCnvs.cnvBed,
      annotatedCnvVcf = annotateCnvs.cnvVcf,
      segments = runControlFreec.segments,
      tumorCpn = tumorCpn.cpn,
      normalCpn = normalCpn.cpn,
      ploidy = runControlFreec.estimatedPloidy,
      tumorPileup = tumorConcatPileups.mergedFile,
      highSensitivityMode = highSensitivityMode,
      gainCutoff = annotateCnvs.gainCutoffUsed,
      lossCutoff = annotateCnvs.lossCutoffUsed,
      genomeModule = genome.genomeModule,
      refFai = genome.refFai,
      image = cnvImage
  }

  output {
    File cnvVcf = annotateCnvs.cnvVcf
    File cnvVcfIndex = annotateCnvs.cnvVcfIndex
    File cnvBed = annotateCnvs.cnvBed
    File filteredCnvBed = filterCnvs.filteredCnvBed
    File segments = runControlFreec.segments
    File freecInfo = runControlFreec.info
    File ratioBedGraph = runControlFreec.ratioBedGraph
    File coveragePlot = filterCnvs.coveragePlot
    File dupDelPlot = filterCnvs.dupDelPlot
    File copyNumberPlot = filterCnvs.copyNumberPlot
    File neutralAfPlot = filterCnvs.neutralAfPlot
    File neutralAfBed = filterCnvs.neutralAfBed
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "Somatic copy-number calling for Ultima Genomics tumor/normal data with [Control-FREEC](https://boevalab.inf.ethz.ch/FREEC/). Per-base coverage of the tumor and the matched normal is collected per chromosome and averaged in 1 kb windows, and pileups at common SNP sites (gnomAD allele frequency above 0.35) give the B-allele frequencies. Control-FREEC segments the tumor/normal coverage ratio and calls gains and losses with allelic support; calls are then annotated, filtered by length and by overlap with Ultima's CNV low-confidence regions, and plotted. Adapted from the Ultima Genomics [controlFREEC_pipeline](https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/controlFREEC_pipeline) workflow (v1.37.1); coverage and pileups are always computed from the input alignments.\n\n`highSensitivityMode` is for low-purity tumors: segments are reported by their fold change (gain above 1.03, loss below 0.97) with ploidy 4 and 10% normal contamination assumed. `chromosomes` restricts the analysis, for a test on part of the genome.\n\nEach task runs its tools from the upstream container images through Apptainer, so the Cromwell backend needs no container support.\n\n## Resources\n\nEverything the workflow reads besides the input alignments comes from three modules, which tasks load through their `modules` runtime attribute:\n\n* `hg38-ultima/v0` sets `HG38_ULTIMA_ROOT`: the hg38 FASTA with its `.fai` and `.dict`.\n* `hg38-ugfreec-resources/1.37.1` sets `HG38_UGFREEC_RESOURCES_ROOT`: the calling regions, 1 kb genome windows, chromosome lengths, common SNP sites, GEM mappability and UG CNV low-confidence regions under `hg38/`, laid out as in the Ultima public bucket.\n* `ugfreec/1.37.1` sets `UGFREEC_ROOT` and loads apptainer: the container images, pulled by digest, under `images/`."
    dependencies: [
      {
        name: "apptainer/1.4.5",
        url: "https://apptainer.org"
      },
      {
        name: "hg38-ultima/v0",
        url: "https://github.com/Ultimagen/healthomics-workflows"
      },
      {
        name: "hg38-ugfreec-resources/1.37.1",
        url: "https://github.com/Ultimagen/healthomics-workflows"
      },
      {
        name: "ultimagenomics/ugbio_freec:1.28.0",
        url: "https://hub.docker.com/r/ultimagenomics/ugbio_freec"
      },
      {
        name: "ultimagenomics/ugbio_cnv:1.33.0",
        url: "https://hub.docker.com/r/ultimagenomics/ugbio_cnv"
      }
    ]
    output_meta: {
      cnvVcf: {
        description: "Copy-number calls with annotations; calls failing the length or low-confidence-region filters are marked in FILTER",
        vidarr_label: "cnvVcf"
      },
      cnvVcfIndex: {
        description: "Tabix index of cnvVcf",
        vidarr_label: "cnvVcfIndex"
      },
      cnvBed: {
        description: "cnvVcf as a BED file, with INFO, FILTER and copy number in the fourth column",
        vidarr_label: "cnvBed"
      },
      filteredCnvBed: {
        description: "cnvBed without the calls that failed a filter",
        vidarr_label: "filteredCnvBed"
      },
      segments: {
        description: "Control-FREEC segmentation of the tumor/normal ratio",
        vidarr_label: "segments"
      },
      freecInfo: {
        description: "Control-FREEC run summary, including the estimated ploidy and normal contamination",
        vidarr_label: "freecInfo"
      },
      ratioBedGraph: {
        description: "Tumor/normal coverage ratio per window, as a BedGraph track",
        vidarr_label: "ratioBedGraph"
      },
      coveragePlot: {
        description: "Normalized coverage along the genome for the tumor and the normal",
        vidarr_label: "coveragePlot"
      },
      dupDelPlot: {
        description: "Gains and losses along the genome",
        vidarr_label: "dupDelPlot"
      },
      copyNumberPlot: {
        description: "Copy number along the genome",
        vidarr_label: "copyNumberPlot"
      },
      neutralAfPlot: {
        description: "Tumor B-allele frequency histogram at SNPs in copy-neutral regions",
        vidarr_label: "neutralAfPlot"
      },
      neutralAfBed: {
        description: "Tumor B-allele frequencies at SNPs in copy-neutral regions",
        vidarr_label: "neutralAfBed"
      }
    }
  }
}

task prepareRegions {
  input {
    Array[String] chromosomes
    String chrLenFile
    String genomeWindows
    String intervalList
    String genomeModule
    Int jobMemory = 2
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  parameter_meta {
    chromosomes: "Chromosomes to keep"
    chrLenFile: "Chromosome lengths, tab-separated name and length"
    genomeWindows: "BED of the genome in equal-sized windows"
    intervalList: "Picard interval list of the calling regions"
    genomeModule: "Modules providing the reference genome and its resource files"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load besides genomeModule"
  }

  command <<<
    set -euo pipefail

    printf '%s\n' ~{sep=" " chromosomes} > requested.txt

    # Keep the genome's own order, and refuse a chromosome the resources do not cover
    awk 'NR == FNR { keep[$1] = 1; next } ($1 in keep)' requested.txt ~{chrLenFile} > chromosomes.genome
    cut -f1 chromosomes.genome > chromosomes.txt
    missing="$(sort requested.txt | comm -23 - <(sort chromosomes.txt))"
    if [[ -n "${missing}" ]]; then
      echo "Not in $(basename ~{chrLenFile}): ${missing}" >&2
      exit 1
    fi

    awk 'NR == FNR { keep[$1] = 1; next } ($1 in keep)' chromosomes.txt ~{genomeWindows} > windows.bed
    # The interval list keeps its whole header: the sequence dictionary describes the reference
    awk 'NR == FNR { keep[$1] = 1; next } /^@/ || ($1 in keep)' chromosomes.txt ~{intervalList} > regions.interval_list
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    Array[String] orderedChromosomes = read_lines("chromosomes.txt")
    File selectedChrLenFile = "chromosomes.genome"
    File selectedWindows = "windows.bed"
    File selectedIntervalList = "regions.interval_list"
  }
}

task scatterIntervals {
  input {
    File intervalList
    Int numShards
    Int breakBandsAtMultiplesOf
    String image
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 1
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    intervalList: "Picard interval list of the regions to call"
    numShards: "Number of shards to split the intervals into"
    breakBandsAtMultiplesOf: "Break intervals at multiples of this many bases"
    image: "Apptainer image providing GATK"
    jobMemory: "Memory allocated for this job (GB); the Java heap is set 3 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    add_bind_root ~{intervalList}
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # With SCATTER_COUNT=1 IntervalListTools writes OUTPUT as a single file rather than
    # a directory of temp_NNNN_of_N subdirectories, so point it at that layout directly.
    mkdir out beds
    output_arg=out
    if [[ ~{numShards} -le 1 ]]; then
      mkdir out/temp_0001_of_1
      output_arg=out/temp_0001_of_1/scattered.interval_list
    fi

    run gatk --java-options "-Xmx~{jobMemory - 3}g" IntervalListTools \
      SCATTER_COUNT=~{numShards} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{breakBandsAtMultiplesOf} \
      INPUT="$(readlink -f ~{intervalList})" \
      OUTPUT="${output_arg}"

    # Number the shards in scatter order so every BED has a unique name
    i=0
    for f in out/*/*.interval_list; do
      i=$((i + 1))
      grep -v '^@' "$f" | awk 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $3 }' > "$(printf 'beds/shard_%04d.bed' "$i")"
    done
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    Array[File] shardBeds = glob("beds/*.bed")
  }
}

task createMpileup {
  input {
    Array[File] crams
    Array[File] cramIndexes
    File shardBed
    String outputPrefix
    Int minMapq
    String genomeModule
    String refFasta
    String snpVcf
    String image
    Int maxDepth = 8000
    Int minBaseQuality = 0
    Int jobMemory = 6
    Int cores = 1
    Int timeout = 4
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    crams: "CRAM or BAM files of one sample"
    cramIndexes: "Indexes of crams, in the same order"
    shardBed: "BED file of the regions in this shard"
    outputPrefix: "Prefix for output file names"
    minMapq: "Minimum mapping quality of reads in the pileup"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA; its .fai and .dict must sit next to it"
    snpVcf: "Bgzipped, indexed VCF of the SNP sites to pile up"
    image: "Apptainer image providing GATK, samtools and bcftools"
    maxDepth: "samtools mpileup maximum depth per file"
    minBaseQuality: "samtools mpileup minimum base quality"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  command <<<
    set -euo pipefail

    crams=(~{sep=" " crams})
    indexes=(~{sep=" " cramIndexes})

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in "${crams[@]}" "${indexes[@]}" ~{shardBed} ~{refFasta} ~{snpVcf}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Cromwell may localize an alignment and its index to different directories; link
    # each pair side by side. The counter prefix keeps files with the same basename apart.
    mkdir alignments tmp
    INPUT_ARGS=""
    for i in "${!crams[@]}"; do
      name="${i}_$(basename "${crams[$i]}")"
      ln -s "${crams[$i]}" "alignments/${name}"
      ln -s "${indexes[$i]}" "alignments/${name}.${indexes[$i]##*.}"
      INPUT_ARGS+=" -I ${PWD}/alignments/${name}"
    done
    REF="$(readlink -f ~{refFasta})"
    BED="$(readlink -f ~{shardBed})"
    SNP="$(readlink -f ~{snpVcf})"
    export INPUT_ARGS REF BED SNP

    # Merge this shard's reads from all files into one BAM, so mpileup reports one sample
    run bash -euo pipefail <<'CONTAINER'
    # shellcheck disable=SC2086
    gatk --java-options "-Xmx2g -Djava.io.tmpdir=${PWD}/tmp" PrintReads ${INPUT_ARGS} \
      -O input.bam -L "${BED}" -R "${REF}"
    samtools index input.bam
    bcftools view -R "${BED}" -Oz -o positions.vcf.gz "${SNP}"
    bcftools index -t positions.vcf.gz
    samtools mpileup -f "${REF}" \
      -d ~{maxDepth} \
      -Q ~{minBaseQuality} \
      -q ~{minMapq} \
      -l positions.vcf.gz \
      input.bam > ~{outputPrefix}.pileup
    rm -f input.bam input.bai input.bam.bai
    CONTAINER
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File pileup = "~{outputPrefix}.pileup"
  }
}

task concatFiles {
  input {
    Array[File] files
    String outputFileName
    Int jobMemory = 2
    Int cores = 1
    Int timeout = 2
    String modules = ""
  }

  parameter_meta {
    files: "Files to concatenate, in order"
    outputFileName: "Name of the concatenated file"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load"
  }

  command <<<
    set -euo pipefail
    cat ~{sep=" " files} > ~{outputFileName}
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File mergedFile = outputFileName
  }
}

task collectCoverage {
  input {
    File cram
    File cramIndex
    String chromosome
    Int minMapq
    String genomeModule
    String refFasta
    String image
    Int jobMemory = 6
    Int cores = 1
    Int timeout = 6
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    cram: "CRAM or BAM file"
    cramIndex: "Index of cram"
    chromosome: "Chromosome to collect coverage over"
    minMapq: "Minimum mapping quality of reads counted"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA, to decode CRAM"
    image: "Apptainer image providing samtools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  String outputFileName = basename(cram) + "." + chromosome + ".bedGraph"

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{cram} ~{cramIndex} ~{refFasta}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # The index is looked up next to the alignment, so pass it by a path that has both
    name="$(basename ~{cram})"
    ln -s ~{cram} "${name}"
    ln -s ~{cramIndex} "${name}.$(basename ~{cramIndex} | awk -F. '{ print $NF }')"
    CRAM="${PWD}/${name}"
    REF="$(readlink -f ~{refFasta})"
    export CRAM REF

    # Per-base depth as a BedGraph, positions included where depth is zero
    run bash -euo pipefail <<'CONTAINER'
    samtools depth -J -Q ~{minMapq} -r ~{chromosome} --reference "${REF}" "${CRAM}" \
      | awk 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $2, $3 }' > ~{outputFileName}
    CONTAINER
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File bedGraph = outputFileName
  }
}

task bedGraphToCpn {
  input {
    Array[File] bedGraphs
    File genomeWindows
    File chrLenFile
    String outputPrefix
    String image
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 4
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    bedGraphs: "Per-base coverage BedGraphs, one per alignment file of the sample"
    genomeWindows: "BED of the genome in equal-sized windows"
    chrLenFile: "Chromosome lengths of the chromosomes to keep"
    outputPrefix: "Prefix for output file names"
    image: "Apptainer image providing bedtools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    : > bedgraphs.txt
    for f in ~{sep=" " bedGraphs} ~{genomeWindows} ~{chrLenFile}; do
      add_bind_root "$f"
    done
    for f in ~{sep=" " bedGraphs}; do
      readlink -f "$f" >> bedgraphs.txt
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    WINDOWS="$(readlink -f ~{genomeWindows})"
    GENOME="$(readlink -f ~{chrLenFile})"
    export WINDOWS GENOME

    # Mean depth per window for each file, windows without coverage set to 0, then the
    # files summed per window. Control-FREEC wants chromosome names without the chr prefix.
    run bash -euo pipefail <<'CONTAINER'
    awk 'BEGIN { OFS = "\t" } { print $1, 0, $2 }' "${GENOME}" > genome.bed
    i=0
    while read -r bedgraph; do
      bedtools intersect -a "${bedgraph}" -b genome.bed -wa > "relevant.${i}.bedgraph"
      bedtools map -g "${GENOME}" -a "${WINDOWS}" -b "relevant.${i}.bedgraph" -c 4 -o mean \
        | awk 'BEGIN { OFS = "\t" } { print $1, $2, $3, ($4 == "." ? 0 : $4) }' \
        | { grep -v chrY || true; } > "mean.${i}.bedgraph"
      rm "relevant.${i}.bedgraph"
      i=$((i + 1))
    done < bedgraphs.txt

    if [[ "${i}" -gt 1 ]]; then
      bedtools unionbedg -i mean.*.bedgraph \
        | awk 'BEGIN { OFS = "\t" } { sum = 0; for (j = 4; j <= NF; j++) sum += $j; print $1, $2, $3, sum }' \
        > union.bedgraph
    else
      cp mean.0.bedgraph union.bedgraph
    fi
    sed 's/^chr//' union.bedgraph > ~{outputPrefix}.cpn
    rm -f mean.*.bedgraph union.bedgraph
    CONTAINER
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File cpn = "~{outputPrefix}.cpn"
  }
}

task runControlFreec {
  input {
    String sampleName
    String tumorCram
    File tumorPileup
    File tumorCpn
    String normalCram
    File normalPileup
    File normalCpn
    File chrLenFile
    String? sex
    Boolean highSensitivityMode
    String genomeModule
    String refFasta
    String refFai
    String snpVcf
    String gemMappability
    String image
    Boolean naiveNormalization = true
    Boolean bedGraphOutput = true
    Boolean contaminationAdjustment = true
    Float? contaminationFraction
    Int? ploidy
    Int? forceGcContentNormalization
    Int window = 1000
    Int degree = 3
    Int jobMemory = 16
    Int cores = 8
    Int timeout = 12
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    sampleName: "Sample name for the Control-FREEC configuration"
    tumorCram: "First tumor alignment file; Control-FREEC names its outputs after it"
    tumorPileup: "Tumor pileup at the SNP sites"
    tumorCpn: "Tumor coverage per window"
    normalCram: "First normal alignment file; Control-FREEC names its outputs after it"
    normalPileup: "Normal pileup at the SNP sites"
    normalCpn: "Normal coverage per window"
    chrLenFile: "Lengths of the chromosomes to call"
    sex: "Sample sex, XX or XY"
    highSensitivityMode: "Assume ploidy 4 and 10% contamination, for low-purity tumors"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA; its .fai must sit next to it"
    refFai: "Reference FASTA index"
    snpVcf: "Bgzipped VCF of the SNP sites"
    gemMappability: "GEM mappability track of the reference"
    image: "Apptainer image providing Control-FREEC and its configuration generator"
    naiveNormalization: "Normalize the tumor by the normal ratio rather than by polynomial fitting"
    bedGraphOutput: "Write the ratios as BedGraph tracks"
    contaminationAdjustment: "Estimate and adjust for normal contamination of the tumor"
    contaminationFraction: "Known normal contamination of the tumor; estimated when unset"
    ploidy: "Expected tumor ploidy; estimated when unset"
    forceGcContentNormalization: "Control-FREEC forceGCcontentNormalization setting: 0, 1 or 2"
    window: "Window size in bases; must match the windows of the coverage files"
    degree: "Degree of the polynomial for GC normalization"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job, also the Control-FREEC thread count"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  String tumorBase = basename(tumorCram)
  String normalBase = basename(normalCram)

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{tumorCram} ~{tumorPileup} ~{tumorCpn} ~{normalCram} ~{normalPileup} ~{normalCpn} \
        ~{chrLenFile} ~{refFasta} ~{refFai} ~{snpVcf} ~{gemMappability}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Control-FREEC writes only the outputs a run produces; create every declared output so
    # a run that has nothing to report still yields the files
    for f in ~{normalBase}_BAF.txt ~{tumorBase}_segments.txt ~{tumorBase}_BAF.txt \
        ~{tumorBase}_CNVs ~{tumorBase}_info.txt ~{tumorBase}_ratio.BedGraph ~{tumorBase}_ratio.txt; do
      touch "$f"
    done

    if ~{highSensitivityMode}; then
      PLOIDY_CONTAMINATION_ARGS="--ploidy 4 --contamination 0.1"
    else
      PLOIDY_CONTAMINATION_ARGS="~{"--ploidy " + ploidy} ~{"--contamination " + contaminationFraction}"
    fi
    TUMOR_CRAM="$(readlink -f ~{tumorCram})"
    NORMAL_CRAM="$(readlink -f ~{normalCram})"
    REF="$(readlink -f ~{refFasta})"
    SNP="$(readlink -f ~{snpVcf})"
    GEM="$(readlink -f ~{gemMappability})"
    CHR_LEN="$(readlink -f ~{chrLenFile})"
    TUMOR_PILEUP="$(readlink -f ~{tumorPileup})"
    NORMAL_PILEUP="$(readlink -f ~{normalPileup})"
    TUMOR_CPN="$(readlink -f ~{tumorCpn})"
    NORMAL_CPN="$(readlink -f ~{normalCpn})"
    export PLOIDY_CONTAMINATION_ARGS TUMOR_CRAM NORMAL_CRAM REF SNP GEM CHR_LEN \
      TUMOR_PILEUP NORMAL_PILEUP TUMOR_CPN NORMAL_CPN

    # Control-FREEC reads the reference one FASTA per chromosome. It takes coverage from the
    # .cpn files and allele counts from the pileups; the alignment paths only name its outputs.
    run bash -euo pipefail <<'CONTAINER'
    mkdir chrFiles_dir
    (cd chrFiles_dir && faidx -x "${REF}" --no-rebuild)

    # shellcheck disable=SC2086
    generate_controlFREEC_config \
      --sample_name ~{sampleName} \
      --BedGraphOutput ~{if bedGraphOutput then "TRUE" else "FALSE"} \
      --NaiveNormalization ~{if naiveNormalization then "TRUE" else "FALSE"} \
      --chrLenFile "${CHR_LEN}" \
      --contaminationAdjustment ~{if contaminationAdjustment then "TRUE" else "FALSE"} \
      --maxThreads ~{cores} \
      --window ~{window} \
      --chrFiles chrFiles_dir \
      --degree ~{degree} \
      ~{"--sex " + sex} \
      --gemMappabilityFile "${GEM}" \
      ~{"--forceGCcontentNormalization " + forceGcContentNormalization} \
      --sample_mateFile "${TUMOR_CRAM}" \
      --sample_mateCopyNumberFile "${TUMOR_CPN}" \
      --sample_miniPileupFile "${TUMOR_PILEUP}" \
      --sample_inputFormat BAM \
      --sample_mateOrientation 0 \
      --control_mateFile "${NORMAL_CRAM}" \
      --control_mateCopyNumberFile "${NORMAL_CPN}" \
      --control_miniPileupFile "${NORMAL_PILEUP}" \
      --control_inputFormat BAM \
      --control_mateOrientation 0 \
      --baf_makePileup "${SNP}" \
      --baf_fastaFile "${REF}" \
      --baf_SNPfile "${SNP}" \
      ${PLOIDY_CONTAMINATION_ARGS}

    /freec -conf ~{sampleName}.config.txt
    rm -rf chrFiles_dir
    CONTAINER

    grep "Output_Ploidy" ~{tumorBase}_info.txt | cut -f2 > ploidy_value.txt
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File segments = "~{tumorBase}_segments.txt"
    File tumorCnvs = "~{tumorBase}_CNVs"
    File info = "~{tumorBase}_info.txt"
    File ratioBedGraph = "~{tumorBase}_ratio.BedGraph"
    Int estimatedPloidy = read_int("ploidy_value.txt")
  }
}

task annotateCnvs {
  input {
    String sampleName
    File cnvCalls
    File segments
    Boolean highSensitivityMode
    String genomeModule
    String refFai
    String cnvLcrBed
    String image
    Float gainCutoff = 1.03
    Float lossCutoff = 0.97
    Int minCnvLength = 10000
    Float intersectionCutoff = 0.5
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 2
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    sampleName: "Sample name in the VCF"
    cnvCalls: "Control-FREEC copy-number calls"
    segments: "Control-FREEC segmentation"
    highSensitivityMode: "Call from segment fold changes rather than from the Control-FREEC calls"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFai: "Reference FASTA index, for the VCF header and the contig naming"
    cnvLcrBed: "BED of CNV low-confidence regions"
    image: "Apptainer image providing the Ultima CNV tools and bcftools"
    gainCutoff: "In highSensitivityMode, minimum segment fold change for a gain"
    lossCutoff: "In highSensitivityMode, maximum segment fold change for a loss"
    minCnvLength: "Calls shorter than this are marked LEN in FILTER"
    intersectionCutoff: "Calls overlapping low-confidence regions by more than this fraction are marked UG-CNV-LCR in FILTER"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  String segmentsBase = basename(segments)
  String cnvsBed = segmentsBase + "_CNVs.bed"
  String outVcf = segmentsBase + "_CNVs.annotate.vcf.gz"

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{cnvCalls} ~{segments} ~{refFai} ~{cnvLcrBed}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # The tools name their outputs after their inputs, so work on links in the task directory
    ln -s ~{segments} ~{segmentsBase}
    CALLS="$(readlink -f ~{cnvCalls})"
    FAI="$(readlink -f ~{refFai})"
    LCR="$(readlink -f ~{cnvLcrBed})"
    export CALLS FAI LCR

    run bash -euo pipefail <<'CONTAINER'
    if ~{highSensitivityMode}; then
      annotate_FREEC_segments \
        --input_segments_file ~{segmentsBase} \
        --gain_cutoff ~{gainCutoff} \
        --loss_cutoff ~{lossCutoff}
    else
      # Neutral calls dropped, contig names brought to the reference's chr convention
      ref_chr=0
      calls_chr=0
      head -1 "${FAI}" | cut -f1 | grep -q '^chr' && ref_chr=1
      head -1 "${CALLS}" | cut -f1 | grep -q '^chr' && calls_chr=1
      if [[ "${ref_chr}" -eq 1 && "${calls_chr}" -eq 0 ]]; then
        sed 's/^/chr/' "${CALLS}" | { grep -v neutral || true; } | cut -f1-4 > ~{cnvsBed}
      elif [[ "${ref_chr}" -eq 0 && "${calls_chr}" -eq 1 ]]; then
        sed 's/^chr//' "${CALLS}" | { grep -v neutral || true; } | cut -f1-4 > ~{cnvsBed}
      else
        { grep -v neutral "${CALLS}" || true; } | cut -f1-4 > ~{cnvsBed}
      fi
    fi

    process_cnvs \
      --sample_name ~{sampleName} \
      --input_bed_file ~{cnvsBed} \
      --intersection_cutoff ~{intersectionCutoff} \
      --min_cnv_length ~{minCnvLength} \
      --out_directory . \
      --fasta_index_file "${FAI}" \
      --cnv_lcr_file "${LCR}"

    if bcftools view -h ~{outVcf} | grep -q '^##FORMAT=<ID=CN,'; then
      bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO;FILTER=%FILTER;CN=[%CN]\n' ~{outVcf} > ~{sampleName}.cnv.bed
    else
      bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO;FILTER=%FILTER\n' ~{outVcf} > ~{sampleName}.cnv.bed
    fi
    CONTAINER
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File cnvVcf = outVcf
    File cnvVcfIndex = outVcf + ".tbi"
    File cnvBed = "~{sampleName}.cnv.bed"
    Float gainCutoffUsed = gainCutoff
    Float lossCutoffUsed = lossCutoff
  }
}

task filterCnvs {
  input {
    String sampleName
    File annotatedCnvBed
    File annotatedCnvVcf
    File segments
    File tumorCpn
    File normalCpn
    Int ploidy
    File tumorPileup
    Boolean highSensitivityMode
    Float gainCutoff
    Float lossCutoff
    String genomeModule
    String refFai
    String image
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 2
    String modules = "ugfreec/1.37.1"
  }

  parameter_meta {
    sampleName: "Sample name for the plots"
    annotatedCnvBed: "Annotated copy-number calls as BED"
    annotatedCnvVcf: "Annotated copy-number calls as VCF"
    segments: "Control-FREEC segmentation"
    tumorCpn: "Tumor coverage per window"
    normalCpn: "Normal coverage per window"
    ploidy: "Ploidy Control-FREEC estimated, the copy number gains and losses are relative to"
    tumorPileup: "Tumor pileup at the SNP sites, for the allele frequency plot"
    highSensitivityMode: "Split gains and losses by fold change rather than by copy number"
    gainCutoff: "In highSensitivityMode, minimum fold change for a gain"
    lossCutoff: "In highSensitivityMode, maximum fold change for a loss"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFai: "Reference FASTA index, for the contig naming"
    image: "Apptainer image providing the Ultima CNV tools and bcftools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugfreec provides apptainer and the container images"
  }

  String filteredBed = basename(segments) + "_CNVs.filter.bed"
  String pileupBase = basename(tumorPileup)

  command <<<
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{annotatedCnvBed} ~{annotatedCnvVcf} ~{tumorCpn} ~{normalCpn} ~{tumorPileup} ~{refFai}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # The allele frequency plot is named after the pileup, so pass it by a link here
    ln -s ~{tumorPileup} ~{pileupBase}
    CNV_BED="$(readlink -f ~{annotatedCnvBed})"
    CNV_VCF="$(readlink -f ~{annotatedCnvVcf})"
    TUMOR_CPN="$(readlink -f ~{tumorCpn})"
    NORMAL_CPN="$(readlink -f ~{normalCpn})"
    FAI="$(readlink -f ~{refFai})"
    # matplotlib needs a writable configuration directory
    MPLCONFIGDIR="${PWD}/matplotlib"
    export CNV_BED CNV_VCF TUMOR_CPN NORMAL_CPN FAI MPLCONFIGDIR

    run bash -euo pipefail <<'CONTAINER'
    { grep -v -E 'FILTER=(UG-CNV-LCR|LEN)' "${CNV_BED}" || true; } > ~{filteredBed}
    bcftools view -f PASS "${CNV_VCF}" | bcftools query -f '%CHROM\t%POS0\t%END\t%CopyNumber\n' > helper.bed

    # Coverage windows as BED, with the chr prefix when the reference uses it
    prefix=""
    head -1 "${FAI}" | cut -f1 | grep -q '^chr' && prefix="chr"
    awk -v p="${prefix}" 'BEGIN { OFS = "\t" } { print p $1, $2, $2 + 999, $NF }' "${NORMAL_CPN}" > normal_coverage.cpn.bed
    awk -v p="${prefix}" 'BEGIN { OFS = "\t" } { print p $1, $2, $2 + 999, $NF }' "${TUMOR_CPN}" > tumor_coverage.cpn.bed

    if ~{highSensitivityMode}; then
      awk '$4 > ~{gainCutoff}' helper.bed > ~{sampleName}.cnvs.filter.DUP.bed
      awk '$4 < ~{lossCutoff}' helper.bed > ~{sampleName}.cnvs.filter.DEL.bed
    else
      awk '$4 > ~{ploidy}' helper.bed > ~{sampleName}.cnvs.filter.DUP.bed
      awk '$4 < ~{ploidy}' helper.bed > ~{sampleName}.cnvs.filter.DEL.bed
    fi

    mkdir -p CNV_figures "${MPLCONFIGDIR}"
    touch CNV_figures/~{sampleName}.CNV.coverage.jpeg \
      CNV_figures/~{sampleName}.dup_del.calls.jpeg \
      CNV_figures/~{sampleName}.CNV.calls.jpeg

    plot_cnv_results \
      --germline_coverage normal_coverage.cpn.bed \
      --tumor_coverage tumor_coverage.cpn.bed \
      --duplication_cnv_calls ~{sampleName}.cnvs.filter.DUP.bed \
      --deletion_cnv_calls ~{sampleName}.cnvs.filter.DEL.bed \
      --sample_name ~{sampleName} \
      --out_directory CNV_figures \
      --neutral_ploidy ~{ploidy}

    plot_FREEC_neutral_AF \
      --mpileup ~{pileupBase} \
      --cnvs_file helper.bed \
      --sample_name ~{sampleName} \
      --out_directory CNV_figures
    CONTAINER

    # Outputs at the top of the task directory, so they are saved without a subdirectory
    mv CNV_figures/* .
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File filteredCnvBed = filteredBed
    File coveragePlot = "~{sampleName}.CNV.coverage.jpeg"
    File dupDelPlot = "~{sampleName}.dup_del.calls.jpeg"
    File copyNumberPlot = "~{sampleName}.CNV.calls.jpeg"
    File neutralAfPlot = "~{pileupBase}.freq.SNP.neutral.hist.jpeg"
    File neutralAfBed = "~{pileupBase}.freq.SNP.neutral.bed"
  }
}
