#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



params.dummy_file = "/home/nehleh/assets/NO_FILE"
dummy_file = file(params.dummy_file)
params.test_file= "/home/nehleh/assets/RMSE_Gubbins.csv"
c_file = file(params.test_file)




frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..5)
nu_sim = Channel.of(0.05) //0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1
recomlen = Channel.of(500) //100,200,300,400,500,1000,2000,3000,4000,5000
recomrate = Channel.of(0.01) //0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1


params.xml = "${PWD}/bin/template/GTR_template.xml"
params.json = "${PWD}/bin/template/GTR_temp_partial.json"

params.genome = 10
params.genomelen = '100000'
params.tMRCA = '0.01'
params.best = false
//params.method = 'pb'
params.hmm_state = '2'
params.nu_hmm = 0.033
params.analyse = false
params.sim_stat = 1 //0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes
params.threshold = 0.5
params.seq = "/home/nehleh/PhyloCode/Result/Results_13092021/num_4/num_4_Wholegenome_4.fasta"
// params.outDir = 'Results'
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
      --best                 true or false (default false)          show the best answer
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)
      --analyse              true or false
  Options for detecting recombination in empirical sequence alignments:
    Mandatory arguments:
      --mode emp
      --seq                  fasta file                             Path to input .fasta file
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)
      --analyse              true or false (default false)

  Output Options:
      --outDir               directory                              Output directory to place final output
  --help                                                            This usage statement

   """.stripIndent()
 }



process MakeClonalTree {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1
     errorStrategy 'ignore'

     input:
         each iteration

     output:
         tuple val(iteration), path('Clonaltree.tree'), emit: Clonaltree
     """
       make_clonaltree.py -n ${params.genome}  -t ${params.tMRCA}
     """
}

process BaciSim {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
         tuple val(iteration), path('Clonaltree')
         each nu_sim
         each recomlen
         each recomrate


     output:
         tuple val(iteration) , path("BaciSimTrees.tree") , emit: BaciSimtrees
         tuple val(iteration) , path('BaciSim_Log.txt') , emit: Recomlog
         tuple val(iteration) , path('BaciSim_Recombination.jpeg'), emit: SimFig
         tuple val(iteration) , path('Sim_nu.csv') , emit: Sim_nu , optional: true
         path('Clonaltree'), emit: Clonaltree
         val iteration , emit: Range
         tuple val(iteration) , path("Recom_stat.csv") , emit: Recomstat
         val nu_sim , emit : nu_sim
         val recomlen , emit : recomlen
         val recomrate , emit: recomrate

     """
       BaciSim.py -cl ${Clonaltree} -n ${params.genome} -g ${params.genomelen} -l ${recomlen} -r ${recomrate}  -nu ${nu_sim} -s ${params.sim_stat}
     """
}


process Seq_gen {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:

        tuple val(iteration), path('BaciSimTrees.tree')
        val f
        val r
        val nu_sim
        val recomlen
        val recomrate

    output:
        path "Wholegenome.fasta" , emit: Wholegenome

    """
     numTrees=\$(wc -l < BaciSimTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomelen} -r$r -f$f -s1 -of BaciSimTrees.tree -p \$numTrees > Wholegenome.fasta
    """
}

process Get_raxml_tree {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Wholegenome
        val iteration
        val nu_sim
        val recomlen
        val recomrate

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
        tuple val(iteration), path('Recomlog')
        path Wholegenome
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'Recom_prob_two.h5'    , emit: recom_prob_two   , optional: true
        path 'PB_Partial_two.xml'   , emit: PB_Partial_two   , optional: true
        path 'PB_Gap_two.xml'       , emit: PB_Gap_two       , optional: true
        path 'PB_Del_two.xml'       , emit: PB_Del_two       , optional: true
        path 'RMSE_PB_two.csv'      , emit :PB_RMSE_two      , optional: true
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_nu_two.txt'        , emit: PB_nu_two        , optional: true
        path 'PB_Partial_eight.xml' , emit: PB_Partial_eight , optional: true
        path 'PB_Gap_eight.xml'     , emit: PB_Gap_eight     , optional: true
        path 'PB_Del_eight.xml'     , emit: PB_Del_eight     , optional: true
        path 'RMSE_PB_eight.csv'    , emit: PB_RMSE_eight    , optional: true
        path 'PB_Recom_eight.jpeg'  , emit: PB_Recom_eight   , optional: true
        path 'PB_Log_eight.txt'     , emit: PB_Log_eight     , optional: true
        path 'PB_nu_eight.txt'      , emit: PB_nu_eight      , optional: true
        path 'PB_Two.catg'          , emit: PB_CATG_two      , optional: true
        path 'PB_two.json'          , emit: PB_JSON_two      , optional: true
        path 'PB_Two_zero.catg'     , emit: PB_CATG_two_zero , optional: true
        path 'PB_two_zero.json'     , emit: PB_JSON_two_zero , optional: true


     """
       phyloHmm_two.py -t ${MyRaxML}  -a ${Wholegenome}  -cl ${Clonaltree} -rl ${Recomlog} -nu ${params.nu_hmm} -st ${params.hmm_state} -sim ${params.simulation} -js ${params.json}

     """
}

process Make_BaciSim_GapDel {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
//     errorStrategy 'ignore'

    input:
        path Clonaltree
        path Wholegenome
        tuple val(iteration), path('Recomlog')
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'Del_alignment.fasta'      , emit: Del_alignment       , optional: true
        path 'Gap_alignment.fasta'      , emit: Gap_alignment       , optional: true
        //path 'BaciSim_partial.json'     , emit: Partial_json        , optional: true
        path 'BaciSim_PSE_999.json'     , emit: BaciSim_PSE_999     , optional: true

        //path 'OriginalSeq.xml' , emit: Original_XML , optional: true
        //path 'BaciSim_Gap.xml' , emit: BaciSim_Gap , optional: true
        //path 'BaciSim_Del.xml' , emit: BaciSim_Del , optional: true
        //path 'BaciSim_partial.xml' , emit: Partial_xml , optional: true
        //path 'BaciSim_partial_certian.xml' , emit: Partial_xml_certian , optional: true
        //path 'BaciSim_Best_Two.catg' , emit: Certain_CATG , optional: true

     """
        BaciSim_GapDel.py -t ${Clonaltree} -a${Wholegenome}  -r${MyRaxML} -l${Recomlog} -x ${params.xml} -j ${params.json}

     """
}


process RaxmlNG_CATG {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
//    errorStrategy 'ignore'

    input:
        path PB_CATG_two
        val iteration
        val nu_sim
        val recomlen
        val recomrate
    output:
        path 'PB_Two.catg.raxml.bestTree', emit: PB_TWO_partialtree

    """
     raxml-ng  --msa ${PB_CATG_two} --msa-format CATG --model GTR+G --prob-msa on
    """
}


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


process PhiloBac_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1
     //errorStrategy 'ignore'

     input:
        path Clonaltree
        tuple val(iteration), path('Recomlog')
        path Wholegenome
        path recom_prob_two
        path physherTree_two
        tuple val(iteration), path ('Recomstat')
        val nu_sim
        val recomlen
        val recomrate


     output:
        path 'RMSE_PB_two.csv'      , emit :PB_RMSE_two      , optional: true
        path 'PB_Recom_two.jpeg'    , emit: PB_Recom_two     , optional: true
        path 'PB_Log_two.txt'       , emit: PB_Log_two       , optional: true
        path 'PB_rcount_two.csv'    , emit: PB_rcount_two    , optional: true
        path 'baci_rcount.csv'      , emit: Baci_rcount      , optional: true
        path 'baci_delta.csv'       , emit: Baci_Delta       , optional: true
        path 'PB_delta_two.csv'     , emit: PB_Delta_two     , optional: true

     """
       PB_result_two.py  -cl ${Clonaltree}  -a ${Wholegenome}  -rl ${Recomlog}  -pb ${physherTree_two} -rp ${recom_prob_two} -st ${params.hmm_state} -sim ${params.simulation} -p ${params.threshold} -rs ${Recomstat}

     """
}




process Get_raxml_tree_gap {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Gap_alignment
        val iteration
        val nu_sim
        val recomlen
        val recomrate
    output:
        path 'RAxML_bestTree.gaptree', emit: GapRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Gap_alignment} -N 10 -n gaptree
    """
}

process Get_raxml_tree_del {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Del_alignment
        val iteration
        val nu_sim
        val recomlen
        val recomrate

    output:
        path 'RAxML_bestTree.deltree', emit: DelRaxML

    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${Del_alignment} -N 10 -n deltree
    """
}


process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
   maxForks 1
//    errorStrategy 'ignore'

     input:
        path Wholegenome
        path MyRaxML
        val iteration
        val nu_sim
        val recomlen
        val recomrate
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
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1
//      errorStrategy 'ignore'

     input:
        path Wholegenome
        path CFML_recom
        path CFMLtree
        tuple val(iteration), path('Recomlog')
        path Clonaltree
        val iteration
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'RMSE_CFML.csv'           , emit : RMSE_CFML  ,optional: true
        path 'CFML_rcount.csv'         , emit: CFML_rcount ,optional: true
        path 'CFML_delta.csv'          , emit: CFML_delta  ,optional: true

     """
       CFML_result.py  -cl ${Clonaltree}  -a ${Wholegenome} -cfl ${CFML_recom}  -cft ${CFMLtree}  -rl ${Recomlog} -sim ${params.simulation}
     """
}


process Run_Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
        path Wholegenome
        val iteration
        val nu_sim
        val recomlen
        val recomrate

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
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
        path Clonaltree
        tuple val(iteration), path('Recomlog')
        path Wholegenome
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat
        val nu_sim
        val recomlen
        val recomrate

    output:
        path 'Gubbins_Recombination.jpeg'    , emit: GubbinsFig
        path 'RMSE_Gubbins.csv'              , emit: Rmse_Gubbins, optional: true
        path 'Gubbinstree_rescale.tree'      , emit: GubbinsRescaletree
        path 'Gubb_rcount.csv'               , emit: Gubb_rcount, optional: true
        path 'Gubbins_delta.csv'             , emit: Gubb_delta, optional: true


    """
     Gubbins_result.py  -cl ${Clonaltree} -a ${Wholegenome}  -rl ${Recomlog}  -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${params.simulation}
    """
}


process Run_Beast {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
         path XML_file
         val prefix
         val iteration
         val nu_sim
         val recomlen
         val recomrate

     output:
         path "${prefix}wholegenome.trees" , emit:  BeastOutput

     """
         beast -prefix ${prefix}  ${XML_file}
     """
}

process Treeannotator {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
         path BeastOutput
         val prefix
         val iteration
         val nu_sim
         val recomlen
         val recomrate

     output:
         path "${prefix}BeastTree.nexus" , emit: BeastNexus

     """
         treeannotator -b 10  ${BeastOutput}  ${prefix}BeastTree.nexus
     """
}


process NexusToNewick {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        path BeastNexus
        val prefix
        val iteration
        val nu_sim
        val recomlen
        val recomrate

     output:
         path "${prefix}BeastTree.newick" , emit : BeastNewick

     """
       NexusToNewick.py -t ${BeastNexus} -o '${prefix}BeastTree.newick'
     """
}


process mergeTreeFiles {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

    input:
         path CFMLtree
         path GubbinsRescaletree
         path physherTree_two
         val iteration
         val nu_sim
         val recomlen
         val recomrate

    output:
         path 'AllOtherTrees.newick' , emit: allOtherTrees

     """
       mergeFiles.py ${CFMLtree} ${GubbinsRescaletree} ${physherTree_two}   > AllOtherTrees.newick
     """
}





process TreeCmp {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
         path Clonaltree
         path allOtherTrees
         val iteration
         val nu_sim
         val recomlen
         val recomrate


     output:
         path 'TreeCmpResult.result' , emit: Comparison

     """
       java -jar /home/nehleh/Documents/0_Research/Software/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${Clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}


process RMSE_summary {

     publishDir "${PWD}/Summary_Results", mode: "copy"
     maxForks 1

     input:
      path  CollectedRMSE_CFML
      path  CollectedRMSE_Gubb
      path  CollectedRMSE_PB_two

     output:
         path 'RMSE_summary.jpeg' , emit: rmse_plot

     """
       rmse_plot_states.py -t rmse_PB_two.csv  -g rmse_Gubbins.csv -c rmse_CFML.csv
     """
}


process RecomCount {

     publishDir "${PWD}/Summary_Results", mode: "copy"
     maxForks 1

     input:
      path  CollectedRcount_Baci
      path  CollectedRcount_PB_two
      path  CollectedRcount_Gubb
      path  CollectedRcount_CFML

     output:
         path 'rcount_scatter.jpeg' , emit: rcount_scatter , optional: true
         path 'rcount_boxplot.jpeg' , emit: rcount_boxplot , optional: true

     """
       rcount_plot.py  -t rcount_PB_two.csv  -b rcount_baci.csv  -g rcount_Gubbins.csv  -c rcount_CFML.csv
     """
}

process Delta_Summary {

     publishDir "${PWD}/Summary_Results", mode: "copy"
     maxForks 1

     input:
      path  CollectedDelta_Baci
      path  CollectedDelta_PB_two
      path  CollectedDelta_Gubb
      path  CollectedDelta_CFML

     output:
         path 'delta_scatter.jpeg' , emit: delta_scatter , optional: true

     """
       delta_plot.py  -t delta_PB_two.csv  -b delta_baci.csv  -c delta_CFML.csv  -g delta_Gubbins.csv
     """
}


process TreeCmp_summary {

     publishDir "${PWD}/Summary_Results", mode: "copy"
     maxForks 1


     input:
        path Comparison


     output:
        path   'TreeCmp_summary.jpeg' , emit: FigTreeCmp

     """
       cmpTree_plot.py -c all_cmpTrees.result
     """
}

workflow Sim {
        take:
            iteration
            frequencies
            rates
            nu_sim
            recomlen
            recomrate
        main:
            MakeClonalTree(iteration)
            BaciSim(MakeClonalTree.out.Clonaltree,nu_sim,recomlen,recomrate)
            Seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates,BaciSim.out.nu_sim,BaciSim.out.recomlen,BaciSim.out.recomrate)

        emit:
            clonaltree = BaciSim.out.Clonaltree
            genome = Seq_gen.out.Wholegenome
            recom_log = BaciSim.out.Recomlog
            iteration = BaciSim.out.Range
            recom_stat = BaciSim.out.Recomstat
            nu_sim = BaciSim.out.nu_sim
            recomlen = BaciSim.out.recomlen
            recomrate= BaciSim.out.recomrate

}


workflow Best {
        take:
            genome
            clonaltree
            raxml_tree
            recom_log
            iteration
            nu_sim
            recomlen
            recomrate
        main:
            Make_BaciSim_GapDel(clonaltree,genome,recom_log,raxml_tree,nu_sim,recomlen,recomrate)
            Get_raxml_tree_gap(Make_BaciSim_GapDel.out.Gap_alignment,iteration,nu_sim,recomlen,recomrate)
            Get_raxml_tree_del(Make_BaciSim_GapDel.out.Del_alignment,iteration,nu_sim,recomlen,recomrate)
        emit:
            GapRaxML = Get_raxml_tree_gap.out.GapRaxML
            DelRaxML = Get_raxml_tree_del.out.DelRaxML
//            Partial_json = Make_BaciSim_GapDel.out.Partial_json
            BaciSim_PSE_999 = Make_BaciSim_GapDel.out.BaciSim_PSE_999
//            original = Make_BaciSim_GapDel.out.Original_XML
//            original_Gap_xml = Make_BaciSim_GapDel.out.BaciSim_Gap
//            original_Del_xml = Make_BaciSim_GapDel.out.BaciSim_Del
//            original_partial_xml = Make_BaciSim_GapDel.out.Partial_xml
//            original_certian_xml = Make_BaciSim_GapDel.out.Partial_xml_certian
//            certian_catg = Make_BaciSim_GapDel.out.Certain_CATG
}


workflow Gubbins {
        take:
            genome
            clonaltree
            recom_log
            iteration
            nu_sim
            recomlen
            recomrate
        main:
            Run_Gubbins(genome,iteration,nu_sim,recomlen,recomrate)
            Gubbins_result(clonaltree,recom_log,genome,Run_Gubbins.out.Gubbinstree,Run_Gubbins.out.GubbinsRecom,Run_Gubbins.out.GubbinsStat,nu_sim,recomlen,recomrate)
        emit:
            RMSE_Gubbins = Gubbins_result.out.Rmse_Gubbins
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree
            Rcount_Gubbins = Gubbins_result.out.Gubb_rcount
            Delta_Gubbins = Gubbins_result.out.Gubb_delta

}

workflow Beast {
        take:
            xml_file
            prefix
            iteration
            nu_sim
            recomlen
            recomrate
        main:
            Run_Beast(xml_file,prefix,iteration,nu_sim,recomlen,recomrate)
            Treeannotator(Run_Beast.out.BeastOutput,prefix,iteration,nu_sim,recomlen,recomrate)
            NexusToNewick(Treeannotator.out.BeastNexus,prefix,iteration,nu_sim,recomlen,recomrate)
        emit:
            BeastNexus = NexusToNewick.out.BeastNewick
}

workflow Beast_cerain {
        take:
            xml_file
            prefix
            iteration
            nu_sim
            recomlen
            recomrate
        main:
            Run_Beast(xml_file,prefix,iteration,nu_sim,recomlen,recomrate)
            Treeannotator(Run_Beast.out.BeastOutput,prefix,iteration,nu_sim,recomlen,recomrate)
            NexusToNewick(Treeannotator.out.BeastNexus,prefix,iteration,nu_sim,recomlen,recomrate)
        emit:
            BeastNexus = NexusToNewick.out.BeastNewick
}



workflow ClonalFrameML {
        take:
            clonaltree
            recom_log
            genome
            raxml_tree
            iteration
            nu_sim
            recomlen
            recomrate

        main:
            CFML(genome,raxml_tree,iteration,nu_sim,recomlen,recomrate)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,recom_log,clonaltree,iteration,nu_sim,recomlen,recomrate)
        emit:
            CFMLtree = CFML.out.CFMLtree
            RMSE_CFML = CFML_result.out.RMSE_CFML
            Rcount_CFML = CFML_result.out.CFML_rcount
            Delta_CFML = CFML_result.out.CFML_delta
}


workflow Physher {
        take:
            json_file
            iteration
            output
            nu_sim
            recomlen
            recomrate

        main:
            Physher_partial(json_file,iteration,nu_sim,recomlen,recomrate)
            Physher_tree(Physher_partial.out.physher_two_txt,iteration,nu_sim,recomlen,recomrate)
        emit:
            physherTree_two = Physher_tree.out.physherTree_two
}



include {   Physher_partial as physher_PSE_9        ; Physher_tree as p_tree_PSE_9 ;
            Physher_partial as physher_PB           ; Physher_tree as p_tree_PB }  from './Physher.nf'



workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.mode == 'sim') {
        println "Detect recombination in simulated data..."
        params.simulation = 1
        Sim(iteration,frequencies,rates,nu_sim,recomlen,recomrate)
        Get_raxml_tree(Sim.out.genome,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)

        if (params.best == true) {
            Best(Sim.out.genome,Sim.out.clonaltree,Get_raxml_tree.out.MyRaxML,Sim.out.recom_log,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
        }
        if (params.method =~ /cfml/) {
            ClonalFrameML(Sim.out.clonaltree,Sim.out.recom_log,Sim.out.genome,Get_raxml_tree.out.MyRaxML,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            CollectedRMSE_CFML = ClonalFrameML.out.RMSE_CFML.collectFile(name:"rmse_CFML.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedRcount_CFML = ClonalFrameML.out.Rcount_CFML.collectFile(name:"rcount_CFML.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedDelta_CFML = ClonalFrameML.out.Delta_CFML.collectFile(name:"delta_CFML.csv",storeDir:"${PWD}/Summary_Results", keepHeader:true , sort: false)

        }
        if (params.method =~ /gub/) {
            Gubbins(Sim.out.genome,Sim.out.clonaltree,Sim.out.recom_log,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            CollectedRMSE_Gubb = Gubbins.out.RMSE_Gubbins.collectFile(name:"rmse_Gubbins.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedRcount_Gubb = Gubbins.out.Rcount_Gubbins.collectFile(name:"rcount_Gubbins.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedDelta_Gubb = Gubbins.out.Delta_Gubbins.collectFile(name:"delta_Gubbins.csv",storeDir:"${PWD}/Summary_Results", keepHeader:true , sort: false)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(Sim.out.clonaltree,Sim.out.recom_log,Sim.out.genome,Get_raxml_tree.out.MyRaxML,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            physher_PB(PhiloBacteria.out.PB_JSON_two,Sim.out.iteration,'physher_PB.txt',Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            p_tree_PB(physher_PB.out.physher_txt,Sim.out.iteration,'physherTree_PB.newick',Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            PhiloBac_result(Sim.out.clonaltree,Sim.out.recom_log,Sim.out.genome,PhiloBacteria.out.recom_prob_two,p_tree_PB.out.physherTree_two,Sim.out.recom_stat,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            CollectedRMSE_PB_two = PhiloBac_result.out.PB_RMSE_two.collectFile(name:"rmse_PB_two.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedRcount_Baci = PhiloBac_result.out.Baci_rcount.collectFile(name:"rcount_baci.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedRcount_PB_two = PhiloBac_result.out.PB_rcount_two.collectFile(name:"rcount_PB_two.csv",storeDir:"${PWD}/Summary_Results", keepHeader:false , sort: false)
            CollectedDelta_Baci = PhiloBac_result.out.Baci_Delta.collectFile(name:"delta_baci.csv",storeDir:"${PWD}/Summary_Results", keepHeader:true , sort: false)
            CollectedDelta_PB_two = PhiloBac_result.out.PB_Delta_two.collectFile(name:"delta_PB_two.csv",storeDir:"${PWD}/Summary_Results", keepHeader:true , sort: false)
        }

         if (params.analyse == true)  {
            RMSE_summary(CollectedRMSE_CFML,CollectedRMSE_Gubb,CollectedRMSE_PB_two)
            RecomCount(CollectedRcount_Baci,CollectedRcount_PB_two,CollectedRcount_Gubb,CollectedRcount_CFML)
            Delta_Summary(CollectedDelta_Baci,CollectedDelta_PB_two,CollectedDelta_Gubb,CollectedDelta_CFML)
            mergeTreeFiles(ClonalFrameML.out.CFMLtree,Gubbins.out.GubbinsRescaletree,p_tree_PB.out.physherTree_two,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            TreeCmp(Sim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            collectedCMP_tree = TreeCmp.out.Comparison.collectFile(name:"all_cmpTrees.result",storeDir:"${PWD}/Summary_Results", keepHeader:true , sort: false)
            TreeCmp_summary(collectedCMP_tree)
         }
    }

     if (params.mode == 'emp') {
        println "Detect recombination in empirical sequence alignments..."
        params.simulation = 0
        recom_log = tuple(1,1,dummy_file)
        clonaltree = dummy_file
        genome = params.seq

        Get_raxml_tree(genome,1,1)

        if (params.method == 'cfml') {
            ClonalFrameML(clonaltree,recom_log,genome,Get_raxml_tree.out.MyRaxML,1,1)
        }
        if (params.method == 'gub') {
            Gubbins(genome,clonaltree,recom_log,1,1)
        }
        if (params.method == 'pb') {

        }

    }

}