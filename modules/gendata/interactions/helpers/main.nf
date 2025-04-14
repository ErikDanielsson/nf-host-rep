process add_interactions_to_phyjson {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple val(genid), path(partial_phyjson_path), path(interactions_csv)

    output:
        tuple val(genid), path("final_phyjson.${genid}.json"), emit: phyjson
    
    script:
    """
    add_interactions_phyjson.py \
        ${partial_phyjson_path} \
        ${interactions_csv} \
        final_phyjson.${genid}.json
    """

    stub:
    """
    touch final_phyjson.${genid}.json
    """

}