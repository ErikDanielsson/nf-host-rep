/*
 * Generate interaction matrices using a TreePPL implementation
 */

process compile_interactions_tppl {
    label 'compile'

    container "${ params.container_treeppl }"

    input: 
        tuple val(param_id), val(seed), path(model_path)
        val(flags) 
    
    output:
        tuple val(param_id), path("sim.${param_id}.out"), emit: sim_bin
    
    script: 
    """
    tpplc $baseDir/bin/simulate.tppl \
        --output sim.${param_id}.out \
        --particles 0 \
        --seed ${seed} \
        ${flags}
    chmod +x sim.${param_id}.out
    """
}

process run_interactions_tppl {
    label 'data'  

    container "${ params.container_treeppl }"
    
    input:
        tuple val(param_id), path(sim_bin), val(genid), path(phyjson_file)

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
            val(param_id), val(genid), path(partial_phyjson_path),
            val(mu), val(beta),
            val(lambda01), val(lambda10), val(lambda12), val(lambda21)
        )
    output:
        tuple val(param_id), val(genid), path("pre_interaction.${param_id}.${genid}.json"), emit: phyjson
    
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
        tuple val(genid), val(param_id), path(int_sim_path), path(host_name_map), path(symbiont_name_map)

    output:
        tuple val(genid), val(param_id), path("interactions.${param_id}.${genid}.csv"), emit: interactions_csv
    
    script:
    """
    interactions_json_to_csv.py \
        ${int_sim_path} \
        ${host_name_map} \
        ${symbiont_name_map} \
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