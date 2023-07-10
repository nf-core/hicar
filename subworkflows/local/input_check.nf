//
// Check input samplesheet and get read channels
//

//
// Include plugins
//
include { fromSamplesheet                          } from 'plugin/nf-validation'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel.fromSamplesheet("input")
        .map{ create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                  // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(ArrayList row) {
    HashMap<String, Object> meta = new HashMap<>(row[0]);
    meta.id           = "${meta.group}_REP${meta.replicate}_T${meta.techniquereplicate}"
    meta.single_end   = false

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row[1]).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row[1]}"
    }
    if (!file(row[2]).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row[2]}"
    }
    array = [ meta, [ file(row[1]), file(row[2]) ] ]

    return array
}
