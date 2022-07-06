#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


frequencies = Channel.value('0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..3)




params.genome = 10
params.genomelen = '3000'
params.recomlen = '500'
params.recomrate = '0.0005'
params.nu_sim = '0.08'
params.nu_hmm = 0.033
params.threshold = 0.9



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
      --analyse              value: true or false
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
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1
     //errorStrategy 'ignore'

     input:
        path Clonaltree
        path Wholegenome
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate
        val iteration

     output:
        path 'PhiloBacter.tree'     ,emit: PB_tree   ,optional: true
        path 'PB_dist.csv'          ,emit: PB_dist   ,optional: true
        path 'PB_WRF_distance.csv'  ,emit: PB_wrf    ,optional: true


     """
       phyloHmm.py -t ${MyRaxML}  -a ${Wholegenome}  -cl ${Clonaltree}  -sim ${params.simulation}

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
        path 'Gubbinstree_rescale.tree'      , emit: GubbinsRescaletree
        path 'Gubb_dist.csv'                 , emit: Gubb_dist            , optional: true
        path 'Gubb_WRF_distance.csv'         , emit: Gubb_wrf             , optional: true


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
            Dist_Gubbins = Gubbins_result.out.Gubb_dist
            WRF_Gubbins = Gubbins_result.out.Gubb_wrf
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
        path 'CFML_dist.csv'           , emit: CFML_dist   ,optional: true
        path 'CFML_WRF_distance.csv'   , emit: CFML_wrf    ,optional: true

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
            Dist_CFML = CFML_result.out.CFML_dist
            WRF_CFML = CFML_result.out.CFML_wrf
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
       cmpTree_plot.py -p dist_PB.csv  -g dist_Gubbins.csv  -m dist_CFML.csv  -prf wrf_PB.csv  -grf  wrf_Gubbins.csv -mrf wrf_CFML.csv
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
            CollectedDist_CFML = ClonalFrameML.out.Dist_CFML.collectFile(name:"dist_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_CFML = ClonalFrameML.out.WRF_CFML.collectFile(name:"wrf_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.method =~ /gub/) {
            Gubbins(Seq_gen.out.Wholegenome,Seq_gen.out.Clonaltree,FastSimBac.out.iteration)
            CollectedDist_Gubb = Gubbins.out.Dist_Gubbins.collectFile(name:"dist_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_Gubb = Gubbins.out.WRF_Gubbins.collectFile(name:"wrf_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(Seq_gen.out.Clonaltree,Seq_gen.out.Wholegenome,Get_raxml_tree.out.MyRaxML,params.nu_sim,params.recomlen,params.recomrate,FastSimBac.out.iteration)
            CollectedDist_PB = PhiloBacteria.out.PB_dist.collectFile(name:"dist_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_PB = PhiloBacteria.out.PB_wrf.collectFile(name:"wrf_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.analyse == true)  {
            TreeCmp_summary(CollectedDist_PB,CollectedDist_Gubb,CollectedDist_CFML,CollectedWRF_PB,CollectedWRF_Gubb,CollectedWRF_CFML)
         }
    }
}