params.datadir = "test_gen_data"

process generate_trees_and_interactions {
    label 'data'

    input:
        val genid
        val nsymbiont
        val nhost

    output:
        tuple val(genid), path("symbiont_tree.${genid}.tre"), emit: symbiont_tree
        tuple val(genid), path("host_tree.${genid}.tre"), emit: host_tree 
        tuple val(genid), path("interactions.${genid}.csv"), emit: interactions_csv 
        tuple val(genid), path("interactions.${genid}.nex"), emit: interactions_nex

    script:
    """
    generate_data.R ${genid} ${nsymbiont} ${nhost}
    """

    stub: 
    """
    touch symbiont_tree.${genid}.tre
    touch host_tree.${genid}.tre
    touch interactions.${genid}.csv
    touch interactions.${genid}.nex
    """
}

process rev_annotate_tree {
    label 'data'

    input:
        val(genid)
        path(input)

    output:
        tuple val(genid), path("${input.getBaseName()}" + ".rev.tre"), emit: rev_tree

    script:
    """
    annotate_tree.Rev --args ${input} --args ${input.baseName}.rev.tre
    """

    stub:
    """
    touch ${input.getBaseName()}.rev.tre
    """
}

process generate_phyjson {
    label 'data'

    input:
        val genid
        path symbiont_tree_file 
        path host_tree_file
        path interactions_csv_file

    output:
        tuple val(genid), path("dirty_host_parasite.${genid}.json"), emit: dirty_phyjson
    
    script:
    """
    transform_data_to_phyjson.R ${symbiont_tree_file} ${host_tree_file} ${interactions_csv_file} dirty_host_parasite.${genid}.json
    """

    stub:
    """
    touch dirty_host_parasite.${genid}.json
    """
}

process clean_phyjson {
    label 'data'

    input:
        val genid
        path dirty_phyjson

    output:
        tuple val(genid), path("host_parasite.${genid}.json"), emit: phyjson

    script:
    """
    clean_phyjson.py ${dirty_phyjson} "host_parasite.${genid}.json"
    """

    stub: 
    """
    touch host_parasite.${genid}.json
    """
}