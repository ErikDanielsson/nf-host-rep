process generate_trees {
    label 'data'

    container "${params.container_r}"

    input:
    val genid
    val nsymbiont
    val nhost

    output:
    tuple val(genid), path("symbiont_tree.${genid}.tre"), emit: symbiont_tree
    tuple val(genid), path("host_tree.${genid}.tre"), emit: host_tree

    script:
    """
    generate_trees.R ${genid} ${nsymbiont} ${nhost}
    """

    stub:
    """
    touch symbiont_tree.${genid}.tre
    touch host_tree.${genid}.tre
    """
}

process rev_annotate_tree {
    label 'data'

    container "${params.container_revbayes}"

    input:
    tuple val(genid), path(input)

    output:
    tuple val(genid), path("${input.getBaseName()}" + ".rev.tre"), emit: rev_tree
    tuple val(genid), path("${input.getBaseName()}" + ".node_index_map.csv"), emit: name_map

    script:
    def labeled_tree_fn = "${input.baseName}.rev.tre"
    def name_map_fn = "${input.baseName}.node_index_map.csv"
    """
    annotate_tree.Rev \
        --args ${input} \
        --args ${labeled_tree_fn}
    node_name_to_index.Rev \
        --args ${labeled_tree_fn} \
        --args ${name_map_fn}
    """

    stub:
    def labeled_tree_fn = "${input.baseName}.rev.tre"
    def name_map_fn = "${input.baseName}.node_index_map.csv"
    """
    touch ${labeled_tree_fn}
    touch ${name_map_fn}
    """
}
/*
process generate_phyjson {
    label 'data'

    container "${ params.container_r }"

    input:
        tuple val(genid), path(symbiont_tree_file), path(host_tree_file), path(interactions_csv_file)

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
*/



process tree_phyjson {
    label 'data'

    container "${params.container_r}"

    input:
    tuple val(genid), path(symbiont_tree_file), path(host_tree_file)

    output:
    tuple val(genid), path("dirty_phyjson.${genid}.json"), emit: dirty_phyjson

    script:
    """
    tree_and_metric_phyjson.R ${symbiont_tree_file} ${host_tree_file} dirty_phyjson.${genid}.json
    """

    stub:
    """
    touch dirty_host_parasite.${genid}.json
    """
}

process clean_phyjson {
    label 'data'

    container "${params.container_python}"

    input:
    tuple val(genid), path(dirty_phyjson)

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
