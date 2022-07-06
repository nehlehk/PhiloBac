#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


frequencies = Channel.value('0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..3)




params.genome = 10
params.genomelen = '3000'
params.recomlen = '500'
params.recomrate = '0.0005'
params.nu_sim = '0.05'
params.outDir = 'SimBacResult'


def helpMessage() {
  log.info"""

  Usage:

  The typical command for running the pipeline is as follows:
    nextflow run main.nf --mode [sim/emp] --other_options

  Options specifying the evolutionary parameters for simulation:
      --mode sim
    Optional arguments:
      --genome               value > 0 (default 10)                 Number of simulated genomes
      --genomelen            value > 1000 (default 100000)          Length of simulated alignment
      --recomlen             value >0 (default 500)                 Length of recombination
      --tMRCA                tMRCA >0 (default 0.01)                TMRCA
      --recomrate            value >0 (default 0.05)                Recombination rate
      --nu_sim               value >0 (default 0.05)                nu value using in simulated data
      --best                 true or false (default)                show the best answer
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)
      --analyse              true or false
  Options for detecting recombination in empirical sequence alignments:
    Mandatory arguments:
      --mode emp
      --seq                  fasta file                             Path to input .fasta file
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)
      --analyse              true or false

  Output Options:
      --outDir               directory                              Output directory to place final output
  --help                                                            This usage statement

   """.stripIndent()
 }


process SimBac {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        each iteration

    output:
        tuple val(iteration) , path("SimBacClonal.tree")   , emit: Clonaltree
        tuple val(iteration) , path("SimBacTrees.tree")    , emit: SimBacTrees
        tuple val(iteration) , path("SimBac_seq.fasta")    , emit: Wholegenome
        tuple val(iteration) , path("SimBac_internal.log") , emit: SimBac_internal
        tuple val(iteration) , path("SimBac_external.log") , emit: SimBac_external
        val iteration                                      , emit: iteration

    """
      SimBac -N ${params.genome}  -B ${params.genomelen} -T 2 -r ${params.recomrate} -e ${params.recomlen} -R 0 -D 0 -m ${params.nu_sim} -M ${params.nu_sim} -c SimBacClonal.tree -l SimBacTrees.tree -b SimBac_internal.log -f SimBac_external.log -o SimBac_seq.fasta
    """
}

process Get_raxml_tree {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
    tuple val(iteration) , path("SimBac_seq")

    output:
        path 'RAxML_bestTree.tree', emit: MyRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${SimBac_seq} -N 10 -n tree
    """
}


process PhiloBacteria {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1
     //errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBacClonal")
        tuple val(iteration) , path("SimBac_seq")
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'PhiloBacter.tree'         ,emit: PB_tree   ,optional: true
        path 'PB_dist.csv'              ,emit: PB_dist   ,optional: true
        path 'PB_WRF_distance.csv'      ,emit: PB_wrf    ,optional: true


     """
       phyloHmm.py -t ${MyRaxML}  -a ${SimBac_seq}  -cl ${SimBacClonal}  -sim ${params.simulation}

     """
}

process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
   maxForks 1
//    errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBac_seq")
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate
    output:
        path "CFML.labelled_tree.newick"    , emit: CFMLtree
        path "CFML.importation_status.txt"  , emit: CFML_recom
        path "CFML.result.txt"              , emit: CFML_result
        path "CFML.em.txt"                  , emit: CFML_em

    """
     ClonalFrameML ${MyRaxML} ${SimBac_seq} CFML >  CFML.result.txt
    """
}

process CFML_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        tuple val(iteration) , path("SimBac_seq")
        tuple val(iteration) , path("SimBacClonal")
        path CFMLtree
        path CFML_recom
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'CFML_dist.csv'           , emit: CFML_dist   ,optional: true
        path 'CFML_WRF_distance.csv'   , emit: CFML_wrf    ,optional: true

     """
       CFML_result.py  -cl ${SimBacClonal}  -a ${SimBac_seq}  -cft ${CFMLtree} -cfl ${CFML_recom}  -sim ${params.simulation}
     """
}

process Run_Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
        tuple val(iteration) , path("SimBac_seq")
        val nu_sim
        val recomlen
        val recomrate

    output:
        path "gubbins.node_labelled.final_tree.tre" , emit: Gubbinstree
        path "gubbins.recombination_predictions.gff" , emit: GubbinsRecom
        path "gubbins.per_branch_statistics.csv" , emit: GubbinsStat
    """
     run_gubbins.py  --model GTRGAMMA  --first-model GTRGAMMA -t raxml  -p gubbins  ${SimBac_seq}
    """
}

process Gubbins_result {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
        tuple val(iteration) , path("SimBacClonal")
        tuple val(iteration) , path("SimBac_seq")
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat
        val nu_sim
        val recomlen
        val recomrate

    output:
        path 'Gubbinstree_rescale.tree'  , emit: GubbinsRescaletree
        path 'Gubb_dist.csv'             , emit: Gubb_dist            , optional: true
        path 'Gubb_WRF_distance.csv'     , emit: Gubb_wrf             , optional: true


    """
     Gubbins_result.py  -cl ${SimBacClonal} -a ${SimBac_seq} -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${params.simulation}
    """
}


process TreeCmp_summary {
     publishDir "${PWD}/Summary_${params.outDir}", mode: "copy"
     maxForks 1


     input:
        path PB_dist
        path Dist_Gubbins
        path Dist_CFML
        path PB_wrf
        path Gubbins_wrf
        path CFML_wrf

     output:
        path   'Dist_summary.jpeg'    , emit: FigTreeDist

     """
       cmpTree_plot.py  -p dist_PB.csv  -g dist_Gubbins.csv  -m dist_CFML.csv  -prf wrf_PB.csv  -grf  wrf_Gubbins.csv -mrf wrf_CFML.csv
     """
}


workflow {

if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.mode == 'sim') {
        println "Detect recombination in simulated data..."
        params.simulation = 2
        SimBac(iteration)
        Get_raxml_tree(SimBac.out.Wholegenome)
        if (params.method =~ /cfml/) {
            CFML(SimBac.out.Wholegenome,Get_raxml_tree.out.MyRaxML,params.nu_sim,params.recomlen,params.recomrate)
            CFML_result(SimBac.out.Wholegenome,SimBac.out.Clonaltree,CFML.out.CFMLtree,CFML.out.CFML_recom,params.nu_sim,params.recomlen,params.recomrate)
            CollectedDist_CFML = CFML_result.out.CFML_dist.collectFile(name:"dist_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_CFML = CFML_result.out.CFML_wrf.collectFile(name:"wrf_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.method =~ /gub/) {
            Run_Gubbins(SimBac.out.Wholegenome,params.nu_sim,params.recomlen,params.recomrate)
            Gubbins_result(SimBac.out.Clonaltree,SimBac.out.Wholegenome,Run_Gubbins.out.Gubbinstree,Run_Gubbins.out.GubbinsRecom,Run_Gubbins.out.GubbinsStat,params.nu_sim,params.recomlen,params.recomrate)
            CollectedDist_Gubb = Gubbins_result.out.Gubb_dist.collectFile(name:"dist_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_Gubb = Gubbins_result.out.Gubb_wrf.collectFile(name:"wrf_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(SimBac.out.Clonaltree,SimBac.out.Wholegenome,Get_raxml_tree.out.MyRaxML,params.nu_sim,params.recomlen,params.recomrate)
            CollectedDist_PB = PhiloBacteria.out.PB_dist.collectFile(name:"dist_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_PB = PhiloBacteria.out.PB_wrf.collectFile(name:"wrf_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.analyse == true)  {
            TreeCmp_summary(CollectedDist_PB,CollectedDist_Gubb,CollectedDist_CFML,CollectedWRF_PB,CollectedWRF_Gubb,CollectedWRF_CFML)
         }
    }
}