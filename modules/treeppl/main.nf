process compile_model {
    label 'compile'
    container "${ params.container_treeppl }"

    input:
        tuple val(compile_id), val(runid), val(model_dir), val(model_key), path(model_path), val(inference_flags)
        path(lib_path)
        val(niter)

    output:
        tuple val(compile_id), path("${model_dir}.${model_key}.${compile_id}.bin"), emit: hostrep_bin
        tuple val(compile_id), path("alignment_*.${compile_id}.html"), optional: true // This is the alignment analysis HTML
    
    script:
    def out_fn = "${model_dir}.${model_key}.${compile_id}.bin"
    def alignment_html = "alignment_${model_dir}.${model_key}.${compile_id}.html"
    def debug_flags = params.debug_treeppl ? (
        "--debug-iterations --debug-alignment-html ${alignment_html}"
        ) : ""
    """
    tpplc ${model_path} \
        --output ${out_fn} \
        --particles ${niter} \
        --seed ${runid} \
        ${debug_flags} \
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
    container "${ params.container_treeppl }"

    input:
        tuple val(compile_id), path(hostrep_bin), val(genid), path(phyjson_file)//, path(lib_path)
        val niter
    
    output:
        tuple val(genid), val(compile_id), path("output.${genid}.${compile_id}.json"), emit: output_json
        tuple val(genid), val(compile_id), path("log.${genid}.${compile_id}.txt"), emit: log
    
    script:
    """
    OCAMLRUNPARAM=b ./${hostrep_bin} ${phyjson_file} ${niter} > output.${genid}.${compile_id}.json 2> log.${genid}.${compile_id}.txt
    """

    stub:
    """
    touch output.${genid}.${compile_id}.json
    echo $niter > log.${genid}.${compile_id}.txt
    """
}
