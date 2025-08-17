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
    --mafft_opts  MAFFT options (default: --globalpair --maxiterate 1000 --thread \$task.cpus)
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

    publishDir "${params.outdir}/busco_dataset", mode: 'copy'

    input:
    val lineage

    output:
    path "${lineage}", emit: lineage_dir

    script:
    """
    busco --download "${lineage}"
    mv "busco_downloads/lineages/${lineage}" "${lineage}"
    """

    stub:
    """
    echo "Stub process for downloading BUSCO dataset: ${lineage}"
    mkdir -p "${lineage}"
    """
}

// // Run BUSCO for each sample in offline mode
// process busco {
//     label 'process_high_memory'
//     tag

//     publishDir

//     input:

//     output:

//     script:
// }

// // Collect per-gene FASTA files from BUSCO outputs
// process collectSeqs {
//     label 'process_low'
//     tag    'collect'

//     publishDir

//     input:

//     output:

//     script:
// }

// // Select shared genes based on completeness fractions
// process selectGenes {
//     label 'process_low'
//     tag    'select'

//     publishDir

//     input:

//     output:

//     script:
// }

// // Align and trim each gene independently
// process alignGenes {
//     label 'process_high'
//     tag

//     publishDir

//     input:

//     output:

//     script:
// }

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

  // Channel setup
  fasta_ch = Channel.fromPath(params.sample, checkIfExists: true)
                    .splitCsv(strip: true, header: true)
                    // .view()

  // Download BUSCO lineage dataset
  lineage_dir = download_busco_dataset(params.lineage)
  lineage_dir.view { "Downloaded BUSCO lineage: ${it}" }
}
