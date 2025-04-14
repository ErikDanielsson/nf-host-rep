/*
 * Generate interaction matrices with the code from the host rep paper
 * NOTE: Since the root repertoires are simulated from a uniform distribution
 * over all repertoires, the resulting reperotires might be erroneous
 */ 
process revbayes_interactions {
    label 'data'
    
    container "${ params.container_revbayes }"

    input:
        tuple val(genid), path(symbiont_tree_file), path(host_tree_file), path(param_file)
    
    output:
        tuple val(genid), path("phylo_interactions.${genid}.rev_csv"), emit: interactions_csv
        tuple val(genid), path("phylo_interactions.${genid}.nex"), emit: interactions_nex
        tuple val(genid), path("settings.${genid}.txt"), emit: settings

    script:
    """
    Simulate.Rev \
        --args ${genid} \
        --args ${param_file} \
        --args ${symbiont_tree_file} \
        --args ${host_tree_file} \
        --args phylo_interactions.${genid}.nex \
        --args phylo_interactions.${genid}.rev_csv \
        --args settings.${genid}.txt
    """

    stub:
    """
    touch phylo_interactions.${genid}.nex
    touch settings.${genid}.txt

    """
}

process clean_rb_csv {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple val(genid), path(rev_csv)
    
    output:
        tuple val(genid), path("phylo_interaction.${genid}.csv"), emit: interactions_csv
    

    script:
    """
    clean_rev_csv.py ${rev_csv} phylo_interaction.${genid}.csv
    """
}