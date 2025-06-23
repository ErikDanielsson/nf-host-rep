process generate_random_interactions {
    label 'data'

    container "${ params.container_r }"

    input:
        val genid
        val nsymbiont
        val nhost
    
    output:
        tuple val(genid), path("interactions.${genid}.csv"), emit: interactions_csv 
        tuple val(genid), path("interactions.${genid}.nex"), emit: interactions_nex

    script:
    """
    interactions.R ${genid} ${nsymbiont} ${nhost}
    """

    stub: 
    """
    touch interactions.${genid}.csv
    touch interactions.${genid}.nex
    """
}