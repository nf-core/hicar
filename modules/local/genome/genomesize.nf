params.options = [:]

process GENOMESIZE {
    tag "$fasta"
    label 'process_low'

    input:
    path fasta

    output:
    val gs      , emit: size

    exec:
    gs = 0.0
    if (params.macs2_gsize) {
        gs = params.macs2_gsize
    } else {
        //genome size remove all N then * 78%
        fasta.withReader{
            String line
            while( line = it.readLine() ){
                if( ! line =~ /^>/ ){
                    l = line.toLowerCase()
                    gs += l.count("a") + l.count("c") + l.count("g") + l.count("t")
                }
            }
        }
        gs *= 0.78
    }
}
