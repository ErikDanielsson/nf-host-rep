process run_hostrep_revbayes {
    label 'rb_sim'

    container "${ params.container_revbayes }"

    input:
        tuple val(runid), val(genid), path(symbiont_tree_file), path(host_tree_file), val(param_id), path(interactions_nex_file)
        val niter
        
    output:
        tuple val(runid), path("out.${param_id}.${genid}.${runid}.logger.log"), emit: clock_log
        tuple val(runid), path("out.${param_id}.${genid}.${runid}.log"), emit: model_log
        tuple val(runid), path("out.${param_id}.${genid}.${runid}.history.txt"), emit: character_summary_log
        tuple val(runid), path("out.${param_id}.${genid}.${runid}.tre"), emit: phy_symbiont_log
   
    script:
    def sampling_period = params.sampling_period ? (
        "${params.sampling_period}"
    ) : "1"
    """
    Infer_test.Rev \
        --args $runid \
        --args $niter \
        --args ${sampling_period} \
        --args $host_tree_file \
        --args $symbiont_tree_file \
        --args $interactions_nex_file \
        --args \$PWD \
        --args out.${param_id}.${genid}.${runid}
    """

    stub:
    """
    touch out.${param_id}.${genid}.${runid}.logger.log
    touch out.${param_id}.${genid}.${runid}.log
    touch out.${param_id}.${genid}.${runid}.history.txt
    touch out.${param_id}.${genid}.${runid}.tre
    """
}