# nf-core/hicar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [11/03/2021]

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
