#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


frequencies = Channel.value('0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..10)




params.genome = 10
params.genomelen = '100000'
params.recomlen = '500'
params.recomrate = '0.0005'
params.nu_sim = '0.02'
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
      --analyse              value: 0,1,2 (default 2)               0:recombination detection, 1: corrected phylogeny, 2: both 0 and 1
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


process FastSimBac {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        each iteration

    output:
        tuple val(iteration) , path("FastSimBacTrees.tree")    , emit: FastSimBacTrees
        val iteration                                          , emit: iteration

    """
      fastSimBac ${params.genome} ${params.genomelen} -r 0 0 -x ${params.recomrate} ${params.recomlen} -T > FastSimBacTrees.tree
    """
}


process Seq_gen {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1
    //errorStrategy 'ignore'

    input:
        tuple val(iteration) , path("FastSimBacTrees.tree")
        val f
        val r

    output:
        path "Wholegenome_${iteration}.fasta" , emit: Wholegenome
        path "Clonaltree.newick"              , emit: Clonaltree

    """
     \$(head -n 1 FastSimBacTrees.tree > Clonaltree.newick )
     numTrees=\$(wc -l < FastSimBacTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomelen} -r$r -f$f -s1 -of FastSimBacTrees.tree -p \$numTrees > Wholegenome_${iteration}.fasta
    """
}


process Get_raxml_tree {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Wholegenome
        val iteration

    output:
        path 'RAxML_bestTree.tree', emit: MyRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Wholegenome} -N 10 -n tree
    """
}


process PhiloBacteria {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
        path Clonaltree
        path Wholegenome
        path MyRaxML
        val iteration
     output:
        path 'Recom_prob_two.h5'    , emit: recom_prob_two   , optional: true
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_nu_two.txt'        , emit: PB_nu_two        , optional: true
        path 'PB_two.json'          , emit: PB_JSON_two      , optional: true




     """
       phyloHmm.py -t ${MyRaxML}  -a ${Wholegenome}  -cl ${Clonaltree}  -nu ${params.nu_hmm} -st ${params.hmm_state} -sim ${params.simulation} -js ${params.json}
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
        path Clonaltree
        path Wholegenome
        path recom_prob_two
        path physherTree_two
        val iteration

     output:
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_Log_two.txt'       , emit: PB_Log_two       , optional: true
        path 'PB_rcount_two.csv'    , emit: PB_rcount_two    , optional: true
        path 'PB_delta_two.csv'     , emit: PB_Delta_two     , optional: true

     """
       PB_result.py  -cl ${Clonaltree}  -a ${Wholegenome}  -pb ${physherTree_two} -rp ${recom_prob_two} -st ${params.hmm_state} -sim ${params.simulation} -p ${params.threshold}

     """
}


process Run_Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1

     input:
        path Wholegenome
        val iteration

    output:
        path "gubbins.node_labelled.final_tree.tre" , emit: Gubbinstree
        path "gubbins.recombination_predictions.gff" , emit: GubbinsRecom
        path "gubbins.per_branch_statistics.csv" , emit: GubbinsStat
    """
     run_gubbins.py  --model GTRGAMMA  --first-model GTRGAMMA -t raxml  -p gubbins   ${Wholegenome}
    """
}

process Gubbins_result {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
    maxForks 1

     input:
        path Clonaltree
        path Wholegenome
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat
        val iteration

    output:
        path 'Gubbins_Recombination.jpeg'    , emit: GubbinsFig
        path 'Gubbinstree_rescale.tree'      , emit: GubbinsRescaletree
        path 'Gubb_rcount.csv'               , emit: Gubb_rcount, optional: true
        path 'Gubbins_delta.csv'             , emit: Gubb_delta, optional: true


    """
     Gubbins_result.py  -cl ${Clonaltree} -a ${Wholegenome} -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${params.simulation}
    """
}


workflow Gubbins {
        take:
            genome
            clonaltree
            iteration
        main:
            Run_Gubbins(genome,iteration)
            Gubbins_result(clonaltree,genome,Run_Gubbins.out.Gubbinstree,Run_Gubbins.out.GubbinsRecom,Run_Gubbins.out.GubbinsStat,iteration)
        emit:
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree
}


process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
   maxForks 1
   errorStrategy 'ignore'

     input:
        path Wholegenome
        path MyRaxML
        val iteration
    output:
        path "CFML.labelled_tree.newick"    , emit: CFMLtree
        path "CFML.importation_status.txt"  , emit: CFML_recom
        path "CFML.result.txt"              , emit: CFML_result
        path "CFML.em.txt"                  , emit: CFML_em
    """
     ClonalFrameML ${MyRaxML} ${Wholegenome} CFML >  CFML.result.txt
    """
}


process CFML_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
        path Wholegenome
        path CFML_recom
        path CFMLtree
        path Clonaltree
        val iteration

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'CFML_rcount.csv'         , emit: CFML_rcount ,optional: true
        path 'CFML_delta.csv'          , emit: CFML_delta  ,optional: true

     """
       CFML_result.py  -cl ${Clonaltree}  -a ${Wholegenome} -cfl ${CFML_recom}  -cft ${CFMLtree}  -sim ${params.simulation}
     """
}


workflow ClonalFrameML {
        take:
            clonaltree
            genome
            raxml_tree
            iteration

        main:
            CFML(genome,raxml_tree,iteration)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,clonaltree,iteration)
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
         path Clonaltree
         path allOtherTrees
         val iteration

     output:
         path 'TreeCmpResult.result' , emit: Comparison

     """
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${Clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}

process TreeCmp_summary {

     publishDir "${PWD}/Summary_FastSimBac", mode: "copy"
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
        FastSimBac(iteration)
        Seq_gen(FastSimBac.out.FastSimBacTrees,frequencies,rates)
        Get_raxml_tree(Seq_gen.out.Wholegenome,FastSimBac.out.iteration)
        if (params.method =~ /cfml/) {
            ClonalFrameML(Seq_gen.out.Clonaltree,Seq_gen.out.Wholegenome,Get_raxml_tree.out.MyRaxML,FastSimBac.out.iteration)
        }
        if (params.method =~ /gub/) {
            Gubbins(Seq_gen.out.Wholegenome,Seq_gen.out.Clonaltree,FastSimBac.out.iteration)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(Seq_gen.out.Clonaltree,Seq_gen.out.Wholegenome,Get_raxml_tree.out.MyRaxML,FastSimBac.out.iteration)
            Physher_partial(PhiloBacteria.out.PB_JSON_two,FastSimBac.out.iteration,'physher_PB.txt')
            Physher_tree(Physher_partial.out.physher_txt,FastSimBac.out.iteration,'physherTree_PB.newick')
            PhiloBac_result(Seq_gen.out.Clonaltree,Seq_gen.out.Wholegenome,PhiloBacteria.out.recom_prob_two,Physher_tree.out.physherTree_two,FastSimBac.out.iteration)
        }
        if (params.analyse == 2)  {
            mergeTreeFiles(ClonalFrameML.out.CFMLtree,Gubbins.out.GubbinsRescaletree,Physher_tree.out.physherTree_two,FastSimBac.out.iteration)
            TreeCmp(Seq_gen.out.Clonaltree,mergeTreeFiles.out.allOtherTrees,FastSimBac.out.iteration)
            collectedCMP_tree = TreeCmp.out.Comparison.collectFile(name:"all_cmpTrees.result",storeDir:"${PWD}/Summary_FastSimBac", keepHeader:true , sort: false)
            TreeCmp_summary(collectedCMP_tree)
         }
    }
}