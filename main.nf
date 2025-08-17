#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  BUSCO Phylogenomics Pipeline
  Nextflow pipeline for BUSCO-based multigene phylogenomics
  Version: 0.1.0
  Author : Akito Shima (ASUQ)
  Email: akito.shima@oist.jp
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters:
    --sample      Path to sample.csv (header: sample,fasta)
    --lineage     BUSCO lineage dataset (e.g. mycoplasmatota_odb12)

  Optional parameters:
    --help        Show this help message
    --outdir      Output directory (default: ./output)
    --fraction    Comma-separated completeness fractions (default: 0.8,0.9,1.0)
    --busco_opts  Extra BUSCO flags (default: "")
    --mafft_opts  MAFFT options (default: --globalpair --maxiterate 1000)
    --trimal_opts trimAl options (default: -automated1)
    --amas_opts   AMAS concat options (default: --in-format fasta --data-type aa --part-format nexus --cores \$task.cpus)
    --iqtree_opts IQ-TREE options (default: -B 1000 -alrt 1000 -m MFP+MERGE -T \$task.cpus)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error "Please provide all required parameters: --sample and --lineage"
}


//-- Processes -----------------------------------------------------------------

// Download BUSCO dataset for offline use
process download_busco_dataset {
    label 'process_single'

    publishDir "${params.outdir}/busco_downloads", mode: 'copy'

    input:
    val lineage

    output:
    path "busco_downloads", emit: lineage_dir

    script:
    """
    busco --download "${lineage}"
    """

    stub:
    """
    echo "Stub process for downloading BUSCO dataset: ${lineage}"
    mkdir -p "busco_downloads/lineages/${lineage}"
    """
}

// Run BUSCO for each sample in offline mode
process busco {
    label 'process_high_memory'
    tag   "${sample}"

    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    path lineage_dir
    val busco_opts

    output:
    path "${sample}", emit: busco_dir

    script:
    """
    busco --in "${fasta}" \
          --lineage_dataset "${params.lineage}" \
          --out "${sample}/busco_output" \
          --mode genome \
          --cpu ${task.cpus} \
          --offline \
          ${busco_opts}
    """

    stub:
    """
    echo "Stub process for BUSCO: ${sample}"
    mkdir -p "${sample}/busco_output"
    echo "Sample: ${sample}" > "${sample}/busco_output/sample_info.txt"
    echo "Lineage: ${lineage_dir}" >> "${sample}/busco_output/sample_info.txt"
    echo "BUSCO run completed for ${sample}"
    """
}

// Collect per-gene FASTA files from BUSCO outputs
process collect_and_select_genes {
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_busco, stageAs: 'busco/*'
    val fractions

    output:
    path 'seqs', emit: seqs_dir
    path "frac*pct_results", emit: frac_results

    script:
    """
    "${projectDir}/bin/busco_multigene_tree.py" collect \
      --input_dir 'busco' --out_dir '.' --cores ${task.cpus}

    "${projectDir}/bin/busco_multigene_tree.py" select \
      --input_dir 'busco' --out_dir '.' \
      --fraction ${fractions} --cores ${task.cpus}
    """

    stub:
    """
    mkdir -p seqs/raw
    """
}

// Align and trim genes
process align_genes {
    label 'process_high'
    tag   { gene }

    publishDir "${params.outdir}/seqs/aligned", mode: 'copy'
    publishDir "${params.outdir}/seqs/trimmed", mode: 'copy'

    input:
    tuple val(gene), path(faa)

    output:
    path "${gene}_aligned.faa", emit: aligned
    path "${gene}_trimmed.faa", emit: trimmed

    script:
    def mafft_opts  = (params.mafft_opts  ?: '')
    def trimal_opts = (params.trimal_opts ?: '')
    """
    mafft ${mafft_opts} --thread ${task.cpus} "${faa}" > "${gene}_aligned.faa"

    trimal ${trimal_opts} \
      -in  "${gene}_aligned.faa" \
      -out "${gene}_trimmed.faa"
    """
}


// // Concatenate alignments and infer phylogenetic trees
// process inferTrees {
//     label 'process_medium'
//     tag

//     publishDir

//     input:

//     output:

//     script:
// }

//-- Workflow ------------------------------------------------------------------
workflow {
  // Parameter parsing
  if (params.help) {
    helpMessage()
    exit 0
  }

  if (params.sample == null || params.lineage == null) {
    missingParametersError()
    exit 1
  }

  // Read fractions (e.g. "0.8,0.9,1.0"), pick smallest for alignment
  frac_list = (params.fraction as String).tokenize(',').collect { it as Double }.sort()
  smallest_pct = (frac_list[0] * 100) as int

  // Channel setup
  fasta_ch = Channel.fromPath(params.sample, checkIfExists: true)
                    .splitCsv(strip: true, header: true)
                    // .view()

  // Download BUSCO lineage dataset
  busco_db = download_busco_dataset(params.lineage)
  // busco_db.view { "Downloaded BUSCO lineage: ${it}" }

  // Run BUSCO for each sample
  busco_results = busco(fasta_ch, busco_db, params.busco_opts)
                      .collect()

  // Collect and select genes from BUSCO results
  busco_genes = collect_and_select_genes(busco_results, params.fraction)

  // Build gene channel for the smallest fraction only -> per-gene alignment input
  min_frac_gene_ch = busco_genes.frac_results
                                .filter { it.name == "frac${smallest_pct}pct_results" }
                                .map { dir -> file("${dir}/frac${smallest_pct}pct_genes.txt") }
                                .view { "Gene list for fraction ${smallest_pct}%: ${it}" }

  // raw dir path from 'seqs' output
  raw_dir_ch = busco_genes.seqs_dir
                          .map { dir -> file("${dir}/raw") }
                          .view { "Raw gene directory: ${it}" }

  // Create (gene, fasta) tuples by combining gene list + raw dir
  gene_ch = min_frac_gene_ch.combine(raw_dir_ch).flatMap { f, rawdir ->
    def genes = new File(f.toString()).readLines().drop(2)
    genes.collect { g -> tuple(g, file("${rawdir}/${g}.faa")) }
  }

  // // Align & trim each gene
  // aligned = align_genes(gene_ch)

  // // Gather all trimmed files to feed each infer task
  // trimmed_all_ch = aligned.trimmed.collect()
}
