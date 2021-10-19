#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.genomeSize = '100000'
params.recom_len = '600'
params.tMRCA = '0.01'
params.nu_sim = '0.05'
params.out = 'Results'


Genome = Channel.value(10)
frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
repeat_range = Channel.value(1..1)
recom_range = Channel.value(1..1)