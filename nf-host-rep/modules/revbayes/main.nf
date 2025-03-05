process run_hostrep_revbayes {
    label 'sim'

    container "${ params.container_revbayes }"

    input:
        tuple val(runid), val(genid), path(symbiont_tree_file), path(host_tree_file), path(interactions_nex_file)
        val niter
        val freq_subsample
        
    output:
        tuple val(runid), path("out.${genid}.${runid}.logger.log"), emit: clock_log
        tuple val(runid), path("out.${genid}.${runid}.log"), emit: model_log
        tuple val(runid), path("out.${genid}.${runid}.history.txt"), emit: character_summary_log
        tuple val(runid), path("out.${genid}.${runid}.tre"), emit: phy_symbiont_log
   
    script:
    """
    Infer_test.Rev \
        --args $runid \
        --args $niter \
        --args $freq_subsample \
        --args $host_tree_file \
        --args $symbiont_tree_file \
        --args $interactions_nex_file \
        --args \$PWD \
        --args out.${genid}.${runid}
    """

    stub:
    """
    touch out.${genid}.${runid}.logger.log
    touch out.${genid}.${runid}.log
    touch out.${genid}.${runid}.history.txt
    touch out.${genid}.${runid}.tre
    """
}