#!/usr/bin/env nextflow

//def sample_ch = []

new File('raw_data/samplesheet.csv').withReader {reader -> 
    //println(reader.readLine())
    //println(reader.readLine())
}

params.input = "raw_data/samplesheet.csv"

workflow{

    sample_ch = Channel.fromPath(params.input)
    .splitCsv(header: true)
    .map{row ->
        def meta = row.subMap('sampleName', 'pairedEnd')
        def r1 = row.Read_1
        def r2 = row.Read_2
        [meta,r1,r2]
        
    }.view()
    
    //sample_ch.view()


}