#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



process Physher_partial {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        path PB_JSON_two
        val iteration
        val output
        val nu_sim
        val recomlen
        val recomrate


     output:
         path output , emit: physher_txt

     """
       physher  ${PB_JSON_two}   >  ${output}
     """
}


process Physher_tree {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        path physher_txt
        val iteration
        val output
        val nu_sim
        val recomlen
        val recomrate

     output:
         path output , emit: physherTree_two

     """
       physher_result.py -t ${physher_txt} -o ${output}
     """
}



