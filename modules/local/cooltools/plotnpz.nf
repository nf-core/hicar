process PLOTNPZ_BY_COOLTOOLS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::cooltools=0.5.1 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' :
        'biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' }"

    input:
    tuple val(meta), path(npz)
    val format

    output:
    tuple val(meta), path("*.${format}")        , emit:plot
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    #!/usr/bin/env python

    import platform
    from textwrap import dedent

    import yaml
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version()
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)

    resolution = ${meta.bin}
    pile = np.load( "$npz" )
    plt.figure(figsize=(6,6))
    norm = LogNorm(    vmin=10**(-1), vmax=10**1)
    im = plt.imshow(
        pile['pileup'],
        cmap='RdBu_r',
        norm = norm
    );

    dim = pile['pileup'].shape
    plt.xticks(np.arange(0,dim[0],5).astype(int),
        (np.arange(0,dim[0],5).astype(int)-dim[0]//2)*resolution//1000,
    )
    plt.yticks(np.arange(0,dim[1],5).astype(int),
        (np.arange(0,dim[1],5).astype(int)-dim[0]//2)*resolution//1000,
    )
    plt.xlabel("offset from pileup center, kb")
    plt.ylabel("offset from pileup center, kb")

    plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7);

    plt.savefig("${prefix}.${format}")
    """
}
