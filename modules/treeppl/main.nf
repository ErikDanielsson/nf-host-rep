process compile_model {
    label 'compile'
    container "${ workflow.containerEngine == 'podman' ? 'docker.io/' : ''}${ params.container_treeppl }"

    input:
        tuple val(compile_id), val(runid), val(model_dir), val(model_key), path(model_path), val(inference_flags)
        path(lib_path)

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
    ls -la ${lib_path} >> ${out_fn}
    chmod +x ${out_fn}
    """
}

process run_hostrep_treeppl {
    label 'sim'
    container "${ workflow.containerEngine == 'podman' ? 'docker.io/' : ''}${ params.container_treeppl }"

    input:
        tuple val(compile_id), path(hostrep_bin), val(genid), path(phyjson_file)//, path(lib_path)
        val niter
    
    output:
        tuple val(genid), val(compile_id), path("output.${genid}.${compile_id}.json"), path("log.${genid}.${compile_id}.txt"), emit: output_json
    
    script:
    """
    ./${hostrep_bin} ${phyjson_file} ${niter} > output.${genid}.${compile_id}.json 2> log.${genid}.${compile_id}.txt
    """

    stub:
    """
    touch output.${genid}.${compile_id}.json
    echo $niter > log.${genid}.${compile_id}.txt
    """
}
