# nf-busco_phylogenomics


[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
<!-- [![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/) -->

## Introduction

**nf-busco_phylogenomics** is a Nextflow pipeline to build phylogenomic tree from BUSCO outputs


## Usage

First, prepare a samplesheet with your input data that looks as follows:

`sample.csv`:

```csv
sample,fasta
species_1,species_1_genome.fasta
```

Each row represents a fasta file.

If your genomes come from nf-annotation, you can build a BUSCO phylogenomics
samplesheet from its master table and `sample.csv`:

```bash
helpers/build_species_representative_sample.py \
   nf_annotation_master.tsv \
   nf_annotation_sample.csv \
   species_representatives.csv
```

The helper selects rows where `Is_Representative` is `yes`, joins them to the
nf-annotation sample CSV by `accession`, and writes `sample,fasta` output for
species representatives. FASTA paths are preserved as written by default and
validated either as written or relative to the nf-annotation sample CSV.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run ASUQ/busco_phylogenomics \
   -profile <docker/singularity/gwdg/oist/viper-cpu/...> \
   --sample sample.csv \
   --lineage <busco_lineage> \
   --outdir <OUTDIR>
```

For representative species generated from nf-annotation, pass the generated
file instead:

```bash
nextflow run ASUQ/busco_phylogenomics \
   -profile <docker/singularity/gwdg/oist/viper-cpu/...> \
   --sample species_representatives.csv \
   --lineage <busco_lineage> \
   --outdir <OUTDIR>
```

`--fraction` accepts exact decimal values such as `0.999`. When the percentage is not an integer, output labels use `_` in place of the decimal point, for example `0.999 -> frac99_9pct_results`.

BUSCO runs with `--tar` by default to reduce inode pressure from high-file-count
BUSCO subdirectories. The collector reads both tarred and untarred BUSCO
sequence output. Use `--busco_tar false` to keep the legacy untarred BUSCO
output layout.

### MPCDF Viper CPU

Launch `viper-cpu` from `viper05i` with Nextflow work and Apptainer cache data
on shared `/ptmp` storage:

```bash
export NXF_APPTAINER_CACHEDIR="/ptmp/$USER/apptainer-cache"

nextflow run ASUQ/busco_phylogenomics \
  -profile viper-cpu \
  -w "/ptmp/$USER/nf-busco_phylogenomics-work" \
  --sample sample.csv \
  --lineage <busco_lineage> \
  --outdir "/ptmp/$USER/nf-busco_phylogenomics-results"
```

Compute tasks use Slurm without a fixed partition, account, or QoS. The
`download_busco_dataset` process remains intact and runs locally with one CPU
and at most 16 GB; the two-CPU local executor permits at most two
internet-dependent tasks at once. Override the cache with
`--apptainer_cache_dir` and the default `apptainer/1.4.3` module with
`--viper_apptainer_module <module/name>` when necessary.

Viper `/ptmp` is not backed up and inactive files are subject to retention
cleanup, so preserve final results elsewhere when required. `/r` is
login-node-only; `/tmp` and generic `$TMPDIR` are unsuitable for Nextflow work
or shared container caches. The reusable `needs_internet` and opt-in
`process_local_scratch` labels are available; the latter uses `$JOB_TMPDIR`.

Run the dependency-locked development test suite with `pixi run test`.

<!--　> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).　-->

<!-- ## Credits

ASUQ/busco_phylogenomics was originally written by ASUQ. -->

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

<!-- ## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md). -->

<!-- ## Citations -->

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use ASUQ/busco_phylogenomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

<!-- An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x). -->
