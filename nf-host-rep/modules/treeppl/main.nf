process compile_model {
    label 'compile'
    container "${ params.container_treeppl }"

    input:
        tuple val(compile_id), val(runid), val(model_dir), val(model_key), path(model_path), val(inference_flags)

    output:
        tuple val(compile_id), path("${model_dir}.${model_key}.${compile_id}.bin"), emit: hostrep_bin
    
    script:
    def out_fn = "${model_dir}.${model_key}.${compile_id}.bin"
    """
    tpplc ${model_path} \
        --output ${out_fn} \
        --seed ${runid} \
        ${inference_flags}
    chmod +x ${out_fn}
    """

    stub:
    def out_fn = "${model_dir}.${model_key}.${compile_id}.bin"
    """
    cat ${model_path}
    echo ${inference_flags} > ${out_fn}
    chmod +x ${out_fn}
    """
}

process run_hostrep_treeppl {
    label 'sim'
    container "${ params.container_treeppl }"

    /*
    The treeppl implementation is light in memory use for most of the
    execution but uses a lot of memory at the end of the execution (upwards of
    10Gb). This is a hacky way of trying many runs at first, and then settling
    for fewer if they collide to much
    */

    input:
        tuple val(compile_id), path(hostrep_bin), val(genid), path(phyjson_file) 
        val niter
    
    output:
        tuple val(genid), val(compile_id), path("output.${genid}.${compile_id}.json"), emit: output_json
    
    script:
    """
    ./${hostrep_bin} ${phyjson_file} ${niter} > output.${genid}.${compile_id}.json
    """

    stub:
    """
    touch output.${genid}.${compile_id}.json
    """
}