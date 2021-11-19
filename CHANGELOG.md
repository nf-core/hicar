# nf-core/hicar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [11/15/2021]

- fix multiple typo.
- change container for dumpreads.nf, mergepeak.nf, mergereads.nf, r1reads.nf and shiftreads.nf
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
- merge the reviewer comments from [#2](https://github.com/nf-core/hicar/pull/2/)
- remove juicer_tools bin file
- remove install_packages.r file
- rename the parepare_circos to circos_prepare
- update documentation
- update multiqc_config.yaml file
- precheck the conditions before doing differential analysis.

### change on [11/02/2021]

- merge the reviewer comments from [#1](https://github.com/nf-core/hicar/pull/1/)
- change the filename from design.csv to test_samplesheet.csv
- change the filename from samplesheet.csv to test_full_samplesheet.csv
- use nf-core repository URL
- update the documentation of README.md
- update the multiqc_config.yaml file format
- remove the regrexp check for replicate in schema_input.json
- update output.md
- update usage.md
- add container for juicer.nf
- update nextflow_schema.json
