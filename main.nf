#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



params.dummy_file = "/home/nehleh/assets/NO_FILE"
dummy_file = file(params.dummy_file)
params.test_file= "/home/nehleh/assets/RMSE_Gubbins.csv"
c_file = file(params.test_file)



// Genome = Channel.value(10)
frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
repeat_range = Channel.value(1..3)
recom_range = Channel.value(1..3)
simulation = Channel.value(1)
params.xml_file = "${PWD}/bin/template/GTR_template.xml"
json_file = "${PWD}/bin/template/GTR_template.json"




params.genome = 10
params.genomelen = '5000'
params.recomlen = '100'
params.recomrate = '0.02'
params.tMRCA = '0.01'
params.nu_sim = '0.05'
params.sim_stat = 0 //0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes
params.sim_fixed = 0 //0 for fixed number and fixed len of recombination and 1 for normal/random way making recombination events.
// params.outDir = 'Results'
params.seq = "/home/nehleh/PhyloCode/Result/Results_13092021/num_4/num_4_Wholegenome_4.fasta"
// params.help = false



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



process MakeClonalTree {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
         each repeat_range

     output:
         tuple val(repeat_range), path('Clonaltree.tree'), emit: Clonaltree
//          val Genome , emit: Genome
//          val repeat_range , emit: Range
     """
       make_clonaltree.py -n ${params.genome}  -t ${params.tMRCA}
     """
}

process BaciSim {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
     maxForks 1

     input:
//          val Genome
         tuple val(repeat_range), path('Clonaltree')
//          val repeat_range
         each recom_range


     output:
         tuple val(repeat_range) ,val(recom_range), path("BaciSimTrees.tree") , emit: BaciSimtrees
         tuple val(repeat_range) ,val(recom_range), path('BaciSim_Log.txt') , emit: Recomlog
         tuple val(repeat_range) ,val(recom_range), path('BaciSim_Recombination.jpeg'), emit: SimFig
         tuple val(repeat_range) ,val(recom_range), path('Sim_nu.csv') , emit: Sim_nu
         path('Clonaltree'), emit: Clonaltree
         val repeat_range , emit: Range
         val recom_range , emit: RecomRange

     """
       BaciSim.py -cl ${Clonaltree} -n ${params.genome} -g ${params.genomelen} -l ${params.recomlen} -r ${params.recomrate}  -nu ${params.nu_sim} -e ${recom_range} -s ${params.sim_stat} -f ${params.sim_fixed}
     """
}


process Seq_gen {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:

        tuple val(repeat_range), val(recom_range),  path('BaciSimTrees.tree')
        val f
        val r

    output:
        path "Wholegenome_${repeat_range}_${recom_range}.fasta" , emit: Wholegenome

    """
     numTrees=\$(wc -l < BaciSimTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomelen} -r$r -f$f -s1 -of BaciSimTrees.tree -p \$numTrees > Wholegenome_${repeat_range}_${recom_range}.fasta
    """
}

process Get_raxml_tree {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Wholegenome
        val repeat_range
        val recom_range

    output:
        path 'RAxML_bestTree.tree', emit: MyRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Wholegenome} -N 10 -n tree
    """
}

process Make_BaciSim_GapDel {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Clonaltree
        path Wholegenome
        tuple val(repeat_range), val(recom_range), path('Recomlog')
        path MyRaxML
        val repeat_range

     output:
        path 'Del_alignment.fasta' , emit: Del_alignment
        path 'Gap_alignment.fasta' , emit: Gap_alignment
        path 'OriginalSeq.xml' , emit: Original_XML
        path 'BaciSim_Gap.xml' , emit: BaciSim_Gap
        path 'BaciSim_Del.xml' , emit: BaciSim_Del
        path 'BaciSim_partial.xml' , emit: Partial_xml
        path 'BaciSim_partial_certian.xml' , emit: Partial_xml_certian
        path 'BaciSim_partial.json' , emit: Partial_json

     """
        BaciSim_GapDel.py -t ${Clonaltree} -a${Wholegenome}  -r${MyRaxML} -l${Recomlog} -x ${params.xml_file} -j ${params.json_file}

     """
}


process Get_raxml_tree_gap {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Gap_alignment
        val repeat_range
        val recom_range
    output:
        path 'RAxML_bestTree.gaptree', emit: GapRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Gap_alignment} -N 10 -n gaptree
    """
}

process Get_raxml_tree_del {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Del_alignment
        val repeat_range
        val recom_range

    output:
        path 'RAxML_bestTree.deltree', emit: DelRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Del_alignment} -N 10 -n deltree
    """
}


process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
   maxForks 1
//    errorStrategy 'ignore'


     input:
        path Wholegenome
        path MyRaxML
        val repeat_range
        val recom_range

    output:
        path "CFML.labelled_tree.newick" , emit: CFMLtree
        path "CFML.importation_status.txt" , emit: CFML_recom
        path "CFML.result.txt" , emit: CFML_result
        path "CFML.em.txt" , emit: CFML_em

    """
     ClonalFrameML ${MyRaxML} ${Wholegenome} CFML >  CFML.result.txt
    """
}



process CFML_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
     maxForks 1
//      errorStrategy 'ignore'

     input:
        path Wholegenome
        path CFML_recom
        path CFMLtree
        tuple val(repeat_range),val(recom_range), path('Recomlog')
        path Clonaltree
        val simulation

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
//         path 'RMSE_CFML.csv' , emit : RMSE_CFML

     """
       CFML_result.py  -cl ${Clonaltree}  -a ${Wholegenome} -cfl ${CFML_recom}  -cft ${CFMLtree}  -rl ${Recomlog} -sim ${simulation}


     """
}


process Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1

     input:
        path Wholegenome
        val repeat_range
        val recom_range

    output:
//         path "gubbins.final_tree.tre" , emit: Gubbinstree
        path "gubbins.node_labelled.final_tree.tre" , emit: Gubbinstree
        path "gubbins.recombination_predictions.gff" , emit: GubbinsRecom
        path "gubbins.per_branch_statistics.csv" , emit: GubbinsStat
    """
     run_gubbins.py  --model GTRGAMMA  --first-model GTRGAMMA -t raxml  -p gubbins   ${Wholegenome}
    """
}

process Gubbins_result {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1

     input:
        path Clonaltree
        tuple val(repeat_range),val(recom_range), path('Recomlog')
        path Wholegenome
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat
        val simulation

    output:
        path 'Gubbins_Recombination.jpeg' , emit: GubbinsFig
//         path 'RMSE_Gubbins.csv' , emit : Rmse_Gubbins
        path 'Gubbinstree_rescale.tree' , emit: GubbinsRescaletree
    """
     Gubbins_result.py  -cl ${Clonaltree} -a ${Wholegenome}  -rl ${Recomlog}  -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${simulation}
    """
}



workflow Sim {

        MakeClonalTree(repeat_range)
        BaciSim(MakeClonalTree.out.Clonaltree,recom_range)
        Seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates)

        emit:
          clonaltree = BaciSim.out.Clonaltree
          genome = Seq_gen.out.Wholegenome
          recom_log = BaciSim.out.Recomlog
          iteration = BaciSim.out.Range
          recomRange = BaciSim.out.RecomRange
}

workflow Rax {
        take:
            genome
            iteration
            recomRange
        main:
            Get_raxml_tree(genome,iteration,recomRange)
        emit:
            raxml_tree = Get_raxml_tree.out.MyRaxML
}

workflow Best {
        take:
            genome
            clonaltree
            raxml_tree
            recom_log
            iteration
            recomRange
        main:
            Make_BaciSim_GapDel(clonaltree,genome,recom_log,raxml_tree,iteration)
            Get_raxml_tree_gap(Make_BaciSim_GapDel.out.Gap_alignment,iteration,recomRange)
            Get_raxml_tree_del(Make_BaciSim_GapDel.out.Del_alignment,iteration,recomRange)
}

workflow Benchmark {
        take:
            genome
            clonaltree
            raxml_tree
            recom_log
            iteration
            recomRange
            simulation
        main:
            CFML(genome,raxml_tree,iteration,recomRange)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,recom_log,clonaltree,simulation)
            Gubbins(genome,iteration,recomRange)
            Gubbins_result(clonaltree,recom_log,genome,Gubbins.out.Gubbinstree,Gubbins.out.GubbinsRecom,Gubbins.out.GubbinsStat,simulation)
        emit:
            CFMLtree = CFML.out.CFMLtree
//             RMSE_CFML = CFML_result.out.RMSE_CFML
//             RMSE_Gubbins = Gubbins_result.out.Rmse_Gubbins
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree

}



workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.mode == 'sim') {
        println "Detect recombination in simulated data..."
        simulation = 1
        Sim()
        Rax(Sim.out.genome,Sim.out.iteration,Sim.out.recomRange)

        if (params.best == true) {
            Best(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange)
        }
        if (params.method == 'cfml') {
            Benchmark(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange,simulation)
        }
        if (params.method == 'gub') {
            Benchmark(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange,simulation)
        }
        if (params.method == 'pb') {
            Benchmark(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange,simulation)
        }
    }

     if (params.mode == 'emp') {
        println "Detect recombination in empirical sequence alignments..."
        simulation = 0
        recom_log = tuple(1,1,dummy_file)
        clonaltree = dummy_file
        genome = params.seq
        Rax(genome,1,1)
        if (params.submode == 'bench') {
            Benchmark(genome,clonaltree,Rax.out.raxml_tree,recom_log,1,1,simulation)
        }
    }

}
