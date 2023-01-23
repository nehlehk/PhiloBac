# PhiloBacter
PhiloBacter is a software package that detects recombination events and reconstructs a phylogenetic tree of bacterial genomes.

In practice, PhiloBacter consists of two tools, the first can be used independently as a recombination detection tool. To simplify the usage of these tools, validate their result, and compare them to the other state-of-the-art methods, a pipeline has been built using Nextflow software.

# Pipeline introduction

The application of this pipeline can be classified into two modes. The first mode is for all users who have sequences of a group of bacteria and are interested to learn about recombination events or want to have an authentic and reliable phylogenetic tree of those bacteria that were not impacted by recombination. These users can choose the analysis mode to use the pipeline. 


### Dependencies
### Instllation
### 1) Analysis mode
In this case, the pipeline consists of two steps. The first step builds an initial tree from the sequence using RAxML. In the second step, suppose the user is only looking for recombination events and their boundaries in the genome or, in other words, looks for the mosaic pattern of the genome. In that case, the pipeline can respond to this request with good accuracy. In addition, the other option is also available if the user is interested in the phylogeny tree. However, the user does not need to set any parameters; if one runs the pipeline using the default settings, all the steps will be done automatically. 

The bonus capability is that this pipeline is not only specific to PhiloBacter. This is possible, if the user also wants to use two other well-known tools in this field, such as [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) and [Gubbins](https://github.com/nickjcroucher/gubbins), and collect the final trees of all three methods. Here is the command to use the pipeline for real datasets.
> ./nextflow main.nf --mode Analysis --seq genome.fasta  --method pb,cfml,gub
### 2) Simulation mode
> Text that is a quote
