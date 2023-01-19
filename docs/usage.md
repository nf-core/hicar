# nf-core/hicar: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/hicar/usage](https://nf-co.re/hicar/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

HiC on Accessible Regulatory DNA (HiCAR) is a feasible tool to study the long-range interactions
of the cis- and trans-regulatory elements (cREs/TREs) especially for the low input
samples and the samples with no available antibodies. To facilitate the analysis of HiCAR
data, the Nextflow platform based nf-core/hicar pipeline can detect the chromatin interaction
in single nucleotide resolution within and among chromosomes.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analzse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
group,replicate,fastq_1,fastq_2
CONTROL,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL,1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL,1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

HiCAR requires paired-end sequencing. Both fastq_1 and fastq_2 must be provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

If md5_1/2 is provided, the pipeline will check the checksums.

A final samplesheet file consisting of two groups of paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
group,replicate,fastq_1,fastq_2,md5_1,md5_2
CONTROL,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,,
CONTROL,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,,
CONTROL,3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,,
TREATMENT,1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz,,
TREATMENT,2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz,,
TREATMENT,3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,,
TREATMENT,3,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz,,
```

| Column      | Description                                                                                                                                                                           |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `group`     | Custom group name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `replicate` | Biological replicates of the samples.                                                                                                                                                 |
| `fastq_1`   | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                            |
| `fastq_2`   | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                            |
| `md5_1`     | Checksum for fastq_1. The checksums of the files will be check to make sure the file is not truncated if provided.                                                                    |
| `md5_2`     | Checksum for fastq_2. The checksums of the files will be check to make sure the file is not truncated if provided.                                                                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Different levels of 3D organization of chromatin

![Schematic representation of the 3D arrangement of chromatin](images/chromatin_organization.png)

### Call A/B compartments and TADs

On a large scale, the arrangement of chromosomes are organised into two compartments labelled A ("active") and B ("inactive").
A/B compartment-associated regions are on the multi-Mb scale and correlate with either open and expression-active chromatin ("A" compartments) or closed and expression inactive chromatin ("B" compartments). A compartments tend to be gene-rich, have high GC-content, contain histone markers for active transcription, and usually displace the interior of the nucleus. The regions in compartment A tend to interact preferentially with A compartment-associated regions than B compartment-associated ones. B compartments, on the other hand, tend to be gene-poor, compact, contain histone markers for gene silencing, and lie on the nuclear periphery.

A topologically associating domain (TAD) is a smaller size genomic region compare to A/B compartments. It is a self-interacting genomic region. Most of the studies indicate TADs regulate gene expression by limiting the enhancer-promoter interaction to each TAD. A number of proteins are known to be associated with TAD formation. The most studied proteins are the protein CTCF and the protein complex cohesin. It has been shown that the TAD boundaries have high levels of CTCF binding and cohesin/lamina shifting edges.

There are multiple available modules to call A/B compartments and TADs.

- For A/B Compartments calling, available tools are 'cooltools', 'hicExplorer', 'Homer', and 'juicer_tools'.
- For TADs calling, available tools are 'cooltools', 'hicexploer', and 'Homer'.

Here is a short introduction about the tools:

- The [`cooltools`](https://github.com/open2c/cooltools) leverages [`cooler`](https://github.com/open2c/cooler/tree/master/cooler) format to enable flexible and reproducible analysis of high-resolution data. `insulation` tool will be used for TADs calling.
- The [`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/) is a set of programs to process, normalize, analyze and visualize Hi-C and cHi-C data. The [`hicFindTADs`](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html) will be used to call TADs.
- The [`Homer`](http://homer.ucsd.edu/homer/interactions2/HiCpca.html) is a software for motif discovery and next-gen sequencing analysis.
- The [`juicer_tools`](https://github.com/aidenlab/juicer) is a platform for analyzing bin sized Hi-C data.

### Call interactions/loops

Chromatin loops (or significant interactions), represent two inter/intra chromosome regions that interact at a high frequency with one another (high reads density in sequence data). Different from HiC, HiCAR data are biased with one ends or both ends in the open chromatin. The [`MAPS`](https://github.com/ijuric/MAPS) are designed to remove the this kind of biases introduced by the ChIP or Tn5-transposition procedure. However, many tools are hesitated to introduce this kind of model-based analysis for interaction analysis since high frequency interactions must happened within the highly opened chromatin regions. Here `nf-core/hicar` provide multiple choices for interactions calling. Available tools are 'MAPS', ['HiC-DC+'](https://doi.org/10.1038/s41467-021-23749-x) and ['peakachu'](https://doi.org/10.1038/s41467-020-17239-9).

In downstream, differential analysis available for called interactions. Available tools are Bioconductor packages such as `edgeR`, and `diffhic`, and [`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/). We borrowed capture Hi-C analysis pipeline from HiCExplorer to do the differential analysis. Different from `edgeR` and `diffhic` pipeline, HiCExplorer pipeline does not require the replicates. A simple differential analysis by set operation are also available.

For annotation, we will use Bioconductor package [`ChIPpeakAnno`](https://bioconductor.org/packages/ChIPpeakAnno/). Please note that, the involved genes are not only distance based annotation. The most of the interaction calling tools are bin-based caller, and the bin size are kilo-base or even more, which make the annotation difficult. For HiCAR data, the R2 reads are Tn5 insertion site of the open chromatin. And the most of the R2 reads will be a anchor of annotion for the gene promoters. We will annotate the interactions by the annotation of called ATAC (R2) peaks located within the interaction regions.

### Call high resolution interactions

The high resolution interaction caller is also available for confident transcription factor enrichment analysis. However, please note that the high resolution interaction caller are not bin-based but peak based analysis, which uses high computational resources.

### Aggregate peak analysis

Aggregate peak analysis (APA) plots the pileup signals detected by high-resolution interaction data. It is a kind of 2 dimension meta-gene analysis. By providing a list of interested genomic coordinates, the pileup signal will present the enrichment between the interactions and the target interested region. Current available tools for APA are [`cooltools`](https://github.com/open2c/cooltools), [`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/) and [`juicer_tools`](https://github.com/aidenlab/juicer)

### Virtual 4C

Circularized Chromosome Conformation Capture (4C) is a powerful technique for studying the interactions of a specific genomic region with the rest of the genome.
Visualize Hi-C data in a virtual 4C (v4c) format can help user to zoom in the interactions within a specific viewpoint. Current available tools for v4c are [`cooltools`](https://github.com/open2c/cooltools), [`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/) and [`trackViewer`](https://doi.org/10.1038/s41592-019-0430-y).

### Available tools

| Tools        | as paramerter | A/B compartments | TADs    | Interactions | Differential analysis | APA     | V4C     |
| :----------- | :------------ | :--------------- | :------ | :----------- | :-------------------- | :------ | :------ |
| cooltools    | cooltools     | &#9745;          | &#9745; | &#9744;      | &#9744;               | &#9745; | &#9745; |
| diffhic      | diffhic       | &#9744;          | &#9744; | &#9744;      | &#9745;               | &#9744; | &#9744; |
| edgeR        | edger         | &#9744;          | &#9744; | &#9744;      | &#9745;               | &#9744; | &#9744; |
| HiC-DC+      | hicdcplus     | &#9744;          | &#9744; | &#9745;      | &#9744;               | &#9744; | &#9744; |
| hicExplorer  | hicexplorer   | &#9745;          | &#9745; | &#9744;      | &#9745;               | &#9745; | &#9745; |
| homer        | homer         | &#9745;          | &#9745; | &#9744;      | &#9744;               | &#9744; | &#9744; |
| MAPS         | maps          | &#9744;          | &#9744; | &#9745;      | &#9744;               | &#9744; | &#9744; |
| juicer_tools | juicebox      | &#9745;          | &#9744; | &#9744;      | &#9744;               | &#9745; | &#9744; |
| peakachu     | peakachu      | &#9744;          | &#9744; | &#9745;      | &#9744;               | &#9744; | &#9744; |
| setOperation | setOperation  | &#9744;          | &#9744; | &#9744;      | &#9745;               | &#9744; | &#9744; |
| trackViewer  | trackviewer   | &#9744;          | &#9744; | &#9744;      | &#9744;               | &#9744; | &#9745; |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/hicar --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/hicar
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/hicar releases page](https://github.com/nf-core/hicar/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/hicar pipeline is failing after multiple re-submissions of the `MAPS_CUT` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[9d/172ca5] NOTE: Process `NFCORE_HICAR:HICAR:MAPS_MULTIENZYME:MAPS_CUT` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_HICAR:HICAR:MAPS_MULTIENZYME:MAPS_CUT (10000)'

Caused by:
    Process `NFCORE_HICAR:HICAR:MAPS_MULTIENZYME:MAPS_CUT (10000)` terminated with an error exit status (137)

Command executed:
    restriction_cut_multipleenzyme.py \
        -f chr22.fa \
        -s GTAC \
        -p 1 \
        -b 10000 \
        -o 10000_CviQI.cut \
        -c 2

Command exit status:
    137

Command output:
    (empty)
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of hiar](https://nf-co.re/hicar/1.0.0/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `MAPS_CUT` process. The quickest way is to search for `process MAPS_CUT` in the [nf-core/hicar Github repo](https://github.com/nf-core/hicar/search?q=process+MAPS_CUT).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/local/maps/cut.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/hicar/blob/master/modules/local/maps/cut.nf#L8).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/hicar/blob/master/conf/base.config#L36-L40) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `MAPS_CUT` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: MAPS_CUT {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

- Error: `The exit status of the task that caused the workflow execution to fail was: null.`

Check the files are readable for the workflow.

- Error: `Session aborted -- Cause: Unable to execute HTTP request: ngi-igenomes.s3.amazonaws.com`

The internet connection reached the limitation. Try to resume the analysis one hour later.

- Error: `PaddingError: Placeholder of length '80' too short in package`

There is no easy answer here. The new `conda` packages should having a longer prefix (255 characters).
The possible solution now is that try to run the pipeline in a shorter folder path, if at all possible.

- Error: `Not a conda environment` or `command not found`

There is something going wrong with the conda environment building.
Just try to remove the conda environment folder and resume the run.

- Error: `unable to load shared object 'work/conda/env-xxxxxx/lib/R/library/rtracklayer/libs/rtracklayer.dylib', dlopen(rtracklayer.dylib, 6) Library not loaded: @rpath/libssl.1.1.dylib`

The openssl installation have issues for `conda`. Try to reinstall it by
`conda activate work/conda/env-xxxxxx && conda install --force-reinstall -y openssl`

- Error: `error Can't locate Statistics/Basic.pm`

The perl-statistics-basic installed in wrong location. Try to reinstall it by
`conda activate work/conda/env-xxxxx && perl -MCPAN -e 'CPAN::install(Statistics::Basic)'`

- `Error in result[[njob]] <- value : attempt to select less than one element in OneIndex`

The error may caused by out of memory (although the error message seems to be unrelated to memory). Try to set `--peak_pair_block` to a smaller number less than 1e9.

### Known issue with Juicer_tools

If you are using [Juicer_tools](https://github.com/aidenlab/juicer/wiki/) with GPU supported, it is not supported by the containers. We are using [Juicer Tools Pre](https://github.com/aidenlab/juicer/wiki/Pre) to create the [hic files](https://doi.org/10.1016/j.cels.2016.07.002) from aligned HiCAR reads.
We recommend having at least 4GB free RAM to generate the hic files.

### Tips for HPC Users

When using an HPC system you should specify the executor matching your system. Check [available executors](https://www.nextflow.io/docs/latest/executor.html) to use the correct executor and parameters.
This instructs Nextflow to submit pipeline tasks as jobs into your HPC workload manager.
Take SLURM workload manager system as an example for the minimal test, this can be done adding the following lines to the `nextflow.config`.

```nextflow
process.executor = 'slurm' // the workload manager name
process.queueSize = 10 // the job queue size
process.clusterOptions = "-J nextFlowHiCAR -p scavenger" // the options, here -p request a specific partition for the resource allocation. It will be different in your cluster.
```

### Useful resources

- [Nextflow pipeline configuration](https://nf-co.re/usage/configuration)

- [Troubleshooting documentation](https://nf-co.re/docs/usage/troubleshooting)
