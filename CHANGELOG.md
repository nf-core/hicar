# nf-core/hicar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [05/03/2022]

- bump version for release.
- update the README.md.

## v1.0dev - [04/25/2022]

- add `totalLinks` parameter for prepare_circos.
- add filters to chromosome names for `hipeak`.
- add parameter `anchor_peaks`.
- Update `MAPS` for new version of `VGAM`.
- add parameter `publish_dir_mode`.
- add circos plot for MAPS loops.
- update the old syntax.
- prepare reads count for hipeak by chromosome.
- fix the bug for R1 peak calling when the peak is out-of-bound ranges.
- fix the bug in calling peaks for `pospoisson regression`
- simplify the fragment peak calling by call peaks only for fragmentation sites
- add unit test for local `subworkflow`
- change R1 reads name to fragment reads in documentation
- use hdf5 to improve the IO efficency
- fix the file path for igv.py
- improve trackviewer.nf
- update gunzip module
- fix a typo in publish_dir of differential analysis.
- fix the empty output of PCA analysis for `multiQC`.
- replace `zcat` by `gunzip -c`
- update documentation of readme
- fix the wrong estimated genome size for unsupported genome.
- fix multiple typo.
- optimize the code for checksums.nf
- add commit number to source code downloaded from github
- update `igv` module
- update `cooler/dump` module
- replace the `macs2` module by `nf-core/macs2` module
- keep bedtools version consistent
- merge rscripts into `nf` files
- change ch_juicer_tool and ch_circos_config into value type channel
- import md5sum from coreutils
- simplify the shell scripts.
- fix a bug for .hic file creation if the resolution is greater than coverage region.
- add documentation for changed source code for MAPS
- remove `paste` command from igv.nf
- fix the bug that annopeaks.r does not creat folder
- change version number output from txt to yml file.
- update citation.md
- change the juicer_tools download on fly
- resolve questions about installation of R packages
- remove juicer_tools bin file
- remove install_packages.r file
- rename the parepare_circos to circos_prepare
- update documentation
- update multiqc_config.yaml file
- precheck the conditions before doing differential analysis.

### change on [11/02/2021]

- add module to covert pair file to bam for visualization
- decrease the memory cost for differential analysis
- add module to create `circos` plot
- add module `igv`
- add module `juicer`
- update QC documentation
- update the memory cost and add ignore `errorStrategy` for `bedtools` sort
- improve memory cost for modules `trackviewer`, `juicer` and `prepare_counts`
- handle multiple errors for `MAPS`
- update the module to prepare the `macs_gsize`
- fix multiple typos in documentation
- change the filename from design.csv to test_samplesheet.csv
- change the filename from samplesheet.csv to test_full_samplesheet.csv
- use nf-core repository URL
- update the multiqc_config.yaml file format
- remove the regrexp check for replicate in schema_input.json
- update output.md
- update usage.md
- add container for juicer.nf
- update nextflow_schema.json
