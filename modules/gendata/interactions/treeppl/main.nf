/*
 * Generate interaction matrices using a TreePPL implementation
 */

process compile_interactions_tppl {
    label 'compile'

    container "${ params.container_treeppl }"

    input: 
        tuple val(genid), path(model_path)
        path(lib_path)
        val(flags) 
    
    output:
        tuple val(genid), path("sim.${genid}.out"), emit: sim_bin
    
    script: 
    """
    tpplc $baseDir/bin/simulate.tppl \
        --output sim.${genid}.out \
        --particles 1 \
        --seed ${genid} \
        ${flags}
    chmod +x sim.${genid}.out
    """
}

process run_interactions_tppl {
    label 'data'  

    container "${ params.container_treeppl }"
    
    input:
        tuple val(genid), path(sim_bin), path(phyjson_file)

    output:
        tuple val(genid), path("simtree_and_interactions.${genid}.json"), emit: interactions_json
        tuple val(genid), path("log.${genid}.txt"), emit: log

    script: 
    """
    OCAMLRUNPARAM=b ./${sim_bin} ${phyjson_file} > simtree_and_interactions.${genid}.json 2> log.${genid}.txt
    """
}

process add_params_phyjson {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple val(genid), path(partial_phyjson_path), path(params_path)

    output:
        tuple val(genid), path("pre_interaction.${genid}.json"), emit: phyjson
    
    script:
    """
    add_params_phyjson.py ${partial_phyjson_path} ${params_path} pre_interaction.${genid}.json
    """

    stub:
    """
    touch pre_interaction.${genid}.json
    """
}

process interactions_json_to_csv {
    label 'data'

    container "${ params.container_python }"

    input:
        tuple val(genid), path(int_sim_path)

    output:
        tuple val(genid), path("interactions.${genid}.csv"), emit: interactions_csv
    
    script:
    """
    interactions_json_to_csv.py \
        ${int_sim_path} \
        interactions.${genid}.csv
    """

    stub:
    """
    touch interactions.${genid}.csv
    """

}

process interactions_csv_to_nex {
    label 'data'

    container "${ params.container_r }"

    input:
        tuple val(genid), path(interactions_csv_path)
    
    output:
        tuple val(genid), path("interactions.${genid}.nex"), emit: interactions_nex
    
    script:
    """
    interactions_csv_to_nex.R ${interactions_csv_path} interactions.${genid}.nex
    """

    stub:
    """
    touch interactions.${genid}.nex
    """
}