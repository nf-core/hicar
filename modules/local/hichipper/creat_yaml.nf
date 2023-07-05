process HICHIPPER_YAML {
    tag "${meta.id}"
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda 'bioconda::multiqc=1.13'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(hicpro), path(peak)
    path resfrags

    output:
    tuple val(meta), path(hicpro), path(peak), path("$prefix.yaml")       , emit:yaml
    path("versions.yml")                                                  , emit:versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    #!/usr/bin/env python

    import platform
    from textwrap import dedent

    import yaml
    import os

    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version()
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)

    hichipper_yaml_content = {
        "peaks": "$peak",
        "resfrags": "$resfrags",
        "hicpro_output": "$hicpro"
    }
    with open("$prefix.yaml", "w") as f:
        yaml.dump(hichipper_yaml_content, f, default_flow_style=False)
    """
}
