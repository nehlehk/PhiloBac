#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



params.genomeSize = '5000'
params.recom_len = '100'
params.tMRCA = '0.01'
params.nu_sim = '0.05'
params.recom_rate = '0.02'
params.sim_stat = 0 //0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes
params.sim_fixed = 1 //0 for fixed number and fixed len of recombination and 1 for normal/random way making recombination events.
params.out = 'Results'
params.xml_file = "${PWD}/bin/template/GTR_template.xml"
params.json_file = "${PWD}/bin/template/GTR_template.json"



Genome = Channel.value(10)
frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
repeat_range = Channel.value(1..1)
recom_range = Channel.value(1..3)




def HelpMessage() {
  log.info"""
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --mode [sim/emp...] --other_options

  Process arguments:
    --out   [str]     Name of output folder. Default 'baseDir/Results'

  1. Generate simulation datasets:
    * Please define evolutionary parameters in main.nf *
    --mode sim

  2. Visualise/summarise simulation outputs (sequence stats, breakpoints):
    --mode sim_v      Summarise output simulation files for sim_bp
    --simdir [str]    Path to dir which contains folder for simulation files (S4_santa). Default 'baseDir/out'

  3. Benchmark recombination detection methods using simulated data:
    --mode div        Move simulated .fasta into subdirs by size. * Use prior to `--mode bm`*

    --mode bm         Detect recombination in simulated datasets and benchmark methods
    --seqn   [int]    Sample size (number of sequences in alignment) to analyse. * Required for `--mode bm` *
    --simdir [str]    Path to dir which contains folder for simulation files (S4_santa). Default 'baseDir/out'

  4. Detect recombination in empirical sequence alignments:
    --mode emp
    --seq [.fasta]    Path to input .fasta file

  5. Calculate classification metrics:
    --mode class       Determine conditions of simulations vs detected recombination
    --rec_path [str]   Path to folder where output for --mode bm is
    --out      [str]   Path to write files to

   """.stripIndent()
 }



process MakeClonalTree {
     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
         val Genome
         each repeat_range

     output:
         tuple val(repeat_range), path('Clonaltree.tree'), emit: Clonaltree
         val Genome , emit: Genome
     """
       make_clonaltree.py -n ${Genome}  -t ${params.tMRCA}
     """
}

process BaciSim {
     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
     maxForks 1

     input:
         val Genome
         tuple val(repeat_range), path('Clonaltree')
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
       BaciSim.py -cl ${Clonaltree} -n ${Genome} -g ${params.genomeSize} -l ${params.recom_len} -r ${params.recom_rate}  -nu ${params.nu_sim} -e ${recom_range} -s ${params.sim_stat} -f ${params.sim_fixed}
     """
}


process Seq_gen {
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
     seq-gen  -mGTR  -l${params.genomeSize} -r$r -f$f -s1 -of BaciSimTrees.tree -p \$numTrees > Wholegenome_${repeat_range}_${recom_range}.fasta
    """
}

process Get_raxml_tree {
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
   publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
     maxForks 1
//      errorStrategy 'ignore'

     input:
        path Wholegenome
        path CFML_recom
        path CFMLtree
        tuple val(repeat_range),val(recom_range), path('Recomlog')
        path Clonaltree

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'RMSE_CFML.csv' , emit : RMSE_CFML

     """
       CMFL_result.py  -cl ${Clonaltree}  -a ${Wholegenome} -cfl ${CFML_recom}  -cft ${CFMLtree}  -rl ${Recomlog}


     """
}


process Gubbins {
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
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
    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_recom_${recom_range}_$filename" }
    maxForks 1

     input:
        path Clonaltree
        tuple val(repeat_range),val(recom_range), path('Recomlog')
        path Wholegenome
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat


    output:
        path 'Gubbins_Recombination.jpeg' , emit: GubbinsFig
        path 'RMSE_Gubbins.csv' , emit : Rmse_Gubbins
        path 'Gubbinstree_rescale.tree' , emit: GubbinsRescaletree
    """
     Gubbins_result.py  -cl ${Clonaltree} -a ${Wholegenome}  -rl ${Recomlog}  -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat}
    """
}



workflow Sim {
        MakeClonalTree(Genome,repeat_range)
        BaciSim(MakeClonalTree.out.Genome,MakeClonalTree.out.Clonaltree,recom_range)
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
        main:
            CFML(genome,raxml_tree,iteration,recomRange)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,recom_log,clonaltree)
            Gubbins(genome,iteration,recomRange)
            Gubbins_result(clonaltree,recom_log,genome,Gubbins.out.Gubbinstree,Gubbins.out.GubbinsRecom,Gubbins.out.GubbinsStat)
        emit:
            CFMLtree = CFML.out.CFMLtree
            RMSE_CFML = CFML_result.out.RMSE_CFML
            RMSE_Gubbins = Gubbins_result.out.Rmse_Gubbins
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree

}



workflow {
    Sim()
    Rax(Sim.out.genome,Sim.out.iteration,Sim.out.recomRange)
    Best(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange)
    Benchmark(Sim.out.genome,Sim.out.clonaltree,Rax.out.raxml_tree,Sim.out.recom_log,Sim.out.iteration,Sim.out.recomRange)
}



// workflow {
//
//
//     if (params.help) {
//         HelpMessage()
//         exit 0
//     }
//
//     if (params.mode == 'sim') {
//         println "Running simulation..."
//
//         MakeClonalTree(Genome,repeat_range)
//         BaciSim(MakeClonalTree.out.Genome,MakeClonalTree.out.Clonaltree,recom_range)
//         Seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates)
//     }
//
//
//     if (params.mode == 'best') {
//         Get_raxml_tree(Seq_gen.out.Wholegenome,BaciSim.out.Range,BaciSim.out.RecomRange)
//         Make_BaciSim_GapDel(BaciSim.out.Clonaltree,Seq_gen.out.Wholegenome,BaciSim.out.Recomlog,Get_raxml_tree.out.MyRaxML,BaciSim.out.Range)
//         Get_raxml_tree_gap(Make_BaciSim_GapDel.out.Gap_alignment,BaciSim.out.Range,BaciSim.out.RecomRange)
//         Get_raxml_tree_del(Make_BaciSim_GapDel.out.Del_alignment,BaciSim.out.Range,BaciSim.out.RecomRange)
//     }
//
//     if (params.mode == 'bench') {
//         println ""
//
//         Get_raxml_tree(Seq_gen.out.Wholegenome,BaciSim.out.Range,BaciSim.out.RecomRange)
//         CFML(Seq_gen.out.Wholegenome,Get_raxml_tree.out.MyRaxML,BaciSim.out.Range,BaciSim.out.RecomRange)
//         CFML_result(Seq_gen.out.Wholegenome,CFML.out.CFML_recom,CFML.out.CFMLtree,BaciSim.out.Recomlog,BaciSim.out.Clonaltree)
//
//         Gubbins(Seq_gen.out.Wholegenome,BaciSim.out.Range,BaciSim.out.RecomRange)
//         Gubbins_result(BaciSim.out.Clonaltree,BaciSim.out.Recomlog,Seq_gen.out.Wholegenome,Gubbins.out.Gubbinstree,Gubbins.out.GubbinsRecom,Gubbins.out.GubbinsStat)
//
//     }
//
//
//     if (params.mode == 'emp') {
//         println "Detect recombination in empirical sequence alignments..."
//     }
//
//     if (params.mode == 'cmp_rec') {
//         println ""
//     }
//
//     if (params.mode == 'cmp_tree') {
//         println ""
//     }
//
//
//
//     }