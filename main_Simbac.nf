#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


frequencies = Channel.value('0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..5)




params.genome = 10
params.genomelen = '100000'
params.recomlen = '500'
params.recomrate = '0.0005'
params.nu_sim = '0.05'
//params.method = 'pb'
params.hmm_state = '2'
params.nu_hmm = 0.033
params.threshold = 0.5
params.json = "${PWD}/bin/template/GTR_temp_partial.json"



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
      --analyse              value: 0,1,2 (default 2)               0:recombination detection, 1: corrected phylogeny, 2: both 0 and 1

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
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBacClonal")
        tuple val(iteration) , path("SimBac_seq")
        path MyRaxML
     output:
        path 'Recom_prob_two.h5'    , emit: recom_prob_two   , optional: true
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_nu_two.txt'        , emit: PB_nu_two        , optional: true
        path 'PB_two.json'          , emit: PB_JSON_two      , optional: true




     """
       phyloHmm_two.py -t ${MyRaxML}  -a ${SimBac_seq}  -cl ${SimBacClonal}  -nu ${params.nu_hmm} -st ${params.hmm_state} -sim ${params.simulation} -js ${params.json}
     """
}

process Physher_partial {

     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1

     input:
        path PB_JSON_two
        val iteration
        val output


     output:
         path output , emit: physher_txt

     """
       physher  ${PB_JSON_two}   >  ${output}
     """
}


process Physher_tree {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename"}
     maxForks 1

     input:
        path physher_txt
        val iteration
        val output

     output:
         path output , emit: physherTree_two

     """
       physher_result.py -t ${physher_txt} -o ${output}
     """
}


process PhiloBac_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBacClonal")
        tuple val(iteration) , path("SimBac_seq")
        path recom_prob_two
        path physherTree_two


     output:
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_Log_two.txt'       , emit: PB_Log_two       , optional: true
        path 'PB_rcount_two.csv'    , emit: PB_rcount_two    , optional: true
        path 'PB_delta_two.csv'     , emit: PB_Delta_two     , optional: true

     """
       PB_result_two.py  -cl ${SimBacClonal}  -a ${SimBac_seq}  -pb ${physherTree_two} -rp ${recom_prob_two} -st ${params.hmm_state} -sim ${params.simulation} -p ${params.threshold}

     """
}


process Run_Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1

     input:
        tuple val(iteration) , path("SimBac_seq")
        val iteration

    output:
        path "gubbins.node_labelled.final_tree.tre" , emit: Gubbinstree
        path "gubbins.recombination_predictions.gff" , emit: GubbinsRecom
        path "gubbins.per_branch_statistics.csv" , emit: GubbinsStat
    """
     run_gubbins.py  --model GTRGAMMA  --first-model GTRGAMMA -t raxml  -p gubbins   ${SimBac_seq}
    """
}

process Gubbins_result {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1

     input:
        tuple val(iteration) , path("SimBacClonal")
        tuple val(iteration) , path("SimBac_seq")
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat

    output:
        path 'Gubbins_Recombination.jpeg'    , emit: GubbinsFig
        path 'Gubbinstree_rescale.tree'      , emit: GubbinsRescaletree
        path 'Gubb_rcount.csv'               , emit: Gubb_rcount, optional: true
        path 'Gubbins_delta.csv'             , emit: Gubb_delta, optional: true


    """
     Gubbins_result.py  -cl ${SimBacClonal} -a ${SimBac_seq} -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${params.simulation}
    """
}


workflow Gubbins {
        take:
            genome
            clonaltree
            iteration
        main:
            Run_Gubbins(genome,iteration)
            Gubbins_result(clonaltree,genome,Run_Gubbins.out.Gubbinstree,Run_Gubbins.out.GubbinsRecom,Run_Gubbins.out.GubbinsStat)
        emit:
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree
}


process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
   maxForks 1
   errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBac_seq")
        path MyRaxML

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
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
        tuple val(iteration) , path("SimBac_seq")
        path CFML_recom
        path CFMLtree
        tuple val(iteration) , path("SimBacClonal")

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'CFML_rcount.csv'         , emit: CFML_rcount ,optional: true
        path 'CFML_delta.csv'          , emit: CFML_delta  ,optional: true

     """
       CFML_result.py  -cl ${SimBacClonal}  -a ${SimBac_seq} -cfl ${CFML_recom}  -cft ${CFMLtree}  -sim ${params.simulation}
     """
}


workflow ClonalFrameML {
        take:
            clonaltree
            genome
            raxml_tree
            iteration

        main:
            CFML(genome,raxml_tree)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,clonaltree)
        emit:
            CFMLtree = CFML.out.CFMLtree
            Rcount_CFML = CFML_result.out.CFML_rcount
            Delta_CFML = CFML_result.out.CFML_delta
}


process mergeTreeFiles {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1

    input:
         path CFMLtree
         path GubbinsRescaletree
         path physherTree_two
         val iteration


    output:
         path 'AllOtherTrees.newick' , emit: allOtherTrees

     """
       mergeFiles.py ${CFMLtree} ${GubbinsRescaletree} ${physherTree_two}   > AllOtherTrees.newick
     """
}



process TreeCmp {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1

     input:

         tuple val(iteration) , path("SimBacClonal")
         path allOtherTrees


     output:
         path 'TreeCmpResult.result' , emit: Comparison

     """
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${SimBacClonal}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}

process TreeCmp_summary {

     publishDir "${PWD}/Summary_SimBac", mode: "copy"
     maxForks 1

     input:
        path Comparison

     output:
        path   'TreeCmp_summary.jpeg' , emit: FigTreeCmp

     """
       cmpTree_plot.py -c all_cmpTrees.result
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
            ClonalFrameML(SimBac.out.Clonaltree,SimBac.out.Wholegenome,Get_raxml_tree.out.MyRaxML,SimBac.out.iteration)
        }
        if (params.method =~ /gub/) {
            Gubbins(SimBac.out.Wholegenome,SimBac.out.Clonaltree,SimBac.out.iteration)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(SimBac.out.Clonaltree,SimBac.out.Wholegenome,Get_raxml_tree.out.MyRaxML)
            Physher_partial(PhiloBacteria.out.PB_JSON_two,SimBac.out.iteration,'physher_PB.txt')
            Physher_tree(Physher_partial.out.physher_txt,SimBac.out.iteration,'physherTree_PB.newick')
            PhiloBac_result(SimBac.out.Clonaltree,SimBac.out.Wholegenome,PhiloBacteria.out.recom_prob_two,Physher_tree.out.physherTree_two)
        }
        if (params.analyse == true)  {
            mergeTreeFiles(ClonalFrameML.out.CFMLtree,Gubbins.out.GubbinsRescaletree,Physher_tree.out.physherTree_two,SimBac.out.iteration)
            TreeCmp(SimBac.out.Clonaltree,mergeTreeFiles.out.allOtherTrees)
            collectedCMP_tree = TreeCmp.out.Comparison.collectFile(name:"all_cmpTrees.result",storeDir:"${PWD}/Summary_SimBac", keepHeader:true , sort: false)
            TreeCmp_summary(collectedCMP_tree)
         }
    }
}