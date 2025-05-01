/*
 * Generate interaction matrices using a TreePPL implementation
 */

process compile_interactions_tppl {
    label 'compile'

    container "${ params.container_treeppl }"

    input: 
        tuple val(genid), path(model_path)
        val(flags) 
    
    output:
        tuple val(genid), path("sim.${genid}.out"), emit: sim_bin
    
    script: 
    """
    tpplc $baseDir/bin/simulate.tppl \
        --output sim.${genid}.out \
        --particles 0 \
        --seed ${genid} \
        ${flags}
    chmod +x sim.${genid}.out
    """
}

process run_interactions_tppl {
    label 'data'  

    container "${ params.container_treeppl }"
    
    input:
        tuple val(genid), path(sim_bin), val(param_id), path(phyjson_file)

    output:
        tuple val(genid), val(param_id), path("simtree_and_interactions.${param_id}.${genid}.json"), emit: interactions_json
        tuple val(genid), val(param_id), path("log.${param_id}.${genid}.txt"), emit: log

    script: 
    """
    OCAMLRUNPARAM=b ./${sim_bin} \
        ${phyjson_file} > simtree_and_interactions.${param_id}.${genid}.json
        2> log.${param_id}.${genid}.txt \
    """
}


process add_params_phyjson {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple(
            val(genid), val(param_id), path(partial_phyjson_path),
            val(mu), val(beta),
            val(lambda01), val(lambda10), val(lambda12), val(lambda21)
        )
    output:
        tuple val(genid), val(param_id), path("pre_interaction.${param_id}.${genid}.json"), emit: phyjson
    
    script:
    """
    add_params_phyjson.py \
        ${partial_phyjson_path} \
        pre_interaction.${param_id}.${genid}.json \
        ${mu} ${beta} \
        ${lambda01} ${lambda10} ${lambda12} ${lambda21}
    """

    stub:
    """
    touch pre_interaction.${param_id}.${genid}.json
    """
}

process interactions_json_to_csv {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple val(genid), val(param_id), path(int_sim_path)

    output:
        tuple val(genid), val(param_id), path("interactions.${param_id}.${genid}.csv"), emit: interactions_csv
    
    script:
    """
    interactions_json_to_csv.py \
        ${int_sim_path} \
        interactions.${param_id}.${genid}.csv
    """

    stub:
    """
    touch interactions.${param_id}.${genid}.csv
    """

}

process interactions_csv_to_nex {
    label 'data'

    container "${ params.container_r }"

    input:
        tuple val(genid), val(param_id), path(interactions_csv_path)
    
    output:
        tuple val(genid), val(param_id), path("interactions.${param_id}.${genid}.nex"), emit: interactions_nex
    
    script:
    """
    interactions_csv_to_nex.R ${interactions_csv_path} interactions.${param_id}.${genid}.nex
    """

    stub:
    """
    touch interactions.${param_id}.${genid}.nex
    """
}