process compile_hostrep_treeppl {
    label 'compile'

    input:
        tuple val(compile_id), val(runid), val(drift_scale), val(gprob)

    output:
        tuple val(compile_id), path("hostrep.${compile_id}.bin"), emit: hostrep_bin
    
    script:
    """
    tpplc $baseDir/models/host_repertoire.tppl \
        -m mcmc-lw-dk \
        --align \
        --cps none \
        --kernel \
        --drift ${drift_scale}\
        --mcmc-lw-gprob ${gprob} \
        --output hostrep.${compile_id}.bin \
        --seed ${runid}
    chmod +x hostrep.${compile_id}.bin
    """

    stub:
    """
    touch hostrep.${compile_id}.bin
    chmod +x hostrep.${compile_id}.bin
    """
}

process compile_model {
    label 'compile'

    // This accepts the compile id,
    // the seed to the algorithm
    // a key for what model to select,
    // and flags specfifying the inference algorithm
    input:
        tuple val(compile_id), val(runid), val(model_key), val(inference_flags)

    output:
        tuple val(compile_id), path("${model_key}.${compile_id}.bin"), emit: hostrep_bin
    
    script:
    def model_fns = [
        original: "host_rep_original.tppl",                          // Original implementation by Mariana Pires Braga
        no_weight: "host_rep_no_weight.tppl",                // Will sample from the proposal distribution (independence model with image restriction)
        rejection_simple: "host_rep_rejection_simple.tppl",          // Removes weight 0.0; resample; (SMC specific) in independence model, replaces with rejection sampling
        rejection_full: "host_rep_rejection_full.tppl",              // Adds rejection sampling with checking if parasite has host at all times
        uniformizaton_simple: "host_rep_uniformization_simple.tppl", // TODO: Uniformization along branches, no adjustment for false repertoires
        uniformizaton_full: "host_rep_uniformization_full.tppl",     // TODO: Uniformization along branches, adjusts for repertoire along branches
    ] // We could just use string interpolation for this, but I think this is less hacky
    def model_fn = model_fns[model_key]
    def out_fn = "${model_key}.${compile_id}.bin"
    """
    tpplc $baseDir/models/${model_fn} \
        --output ${out_fn} \
        --seed ${runid} \
        ${inference_flags}
    chmod +x ${out_fn}
    """

    stub:
    def model_fns = [
        original: "host_rep_original.tppl",                          // Original implementation by Mariana Pires Braga
        no_weight: "host_rep_no_weight.tppl",                // Will sample from the proposal distribution (independence model with image restriction)
        rejection_simple: "host_rep_rejection_simple.tppl",          // Removes weight 0.0; resample; (SMC specific) in independence model, replaces with rejection sampling
        rejection_full: "host_rep_rejection_full.tppl",              // Adds rejection sampling with checking if parasite has host at all times
        uniformizaton_simple: "host_rep_uniformization_simple.tppl", // TODO: Uniformization along branches, no adjustment for false repertoires
        uniformizaton_full: "host_rep_uniformization_full.tppl",     // TODO: Uniformization along branches, adjusts for repertoire along branches
    ] // We could just use string interpolation for this, but I think this is less hacky
    def model_fn = model_fns[model_key]
    def out_fn = "${model_key}.${compile_id}.bin"
    """
    cat $baseDir/models/${model_fn}
    echo ${inference_flags} > ${out_fn}
    chmod +x ${out_fn}
    """
}

process run_hostrep_treeppl {
    label 'sim'
    /*
    The treeppl implementation is light in memory use for most of the
    execution but uses a lot of memory at the end of the execution (upwards of
    10Gb). This is a hacky way of trying many runs at first, and then settling
    for fewer if they collide to much
    */
    memory { 1.5.GB * Math.pow(2, task.attempt - 1) }
    maxRetries 5 
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

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

process time_hostrep_treeppl {
    label 'sim'
    memory { 2.GB * Math.pow(2, task.attempt - 1) }
    maxRetries 5 

    // publishDir "${params.outdir}"

    input:
        tuple val(runid), path(hostrep_bin), val(genid), path(phyjson_file) 
        val niter

    output:
        tuple val(runid), path("time.treeppl.${genid}.${runid}.txt"), emit: time
        tuple val(runid), path("output.${genid}.${runid}.json"), emit: output_json
    
    script:
    """
    { time \
        ./${hostrep_bin} ${phyjson_file} ${niter} \
        1> output.${genid}.${runid}.json \
        2> /dev/null; \
    } 2> "time.treeppl.${genid}.${runid}.txt"
    """


}

process perf_hostrep_treeppl {
    publishDir "${params.outdir}"

    input:
        val runid
        val niter
        path hostrep_bin
        path phyjson_file
    
    output:
        path "time.treeppl.${runid}.txt"
        path "output.${runid}.json" 
    
    script:
    """
    { perf \
        ./${hostrep_bin} ${phyjson_file} ${niter} \
        1> output.${runid}.json \
        2> /dev/null; \
    } 2> "time.treeppl.${runid}.txt"
    """
}