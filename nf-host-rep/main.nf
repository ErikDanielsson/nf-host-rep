nextflow.enable.dsl=2

include {
    generate_trees_and_interactions;
    rev_annotate_tree;
    generate_phyjson;
    clean_phyjson
} from "./modules/gendata"

include {
    compile_model;
    run_hostrep_treeppl;
} from "./modules/treeppl"

include {
    run_hostrep_revbayes
} from "./modules/revbayes"

workflow {
    // Define the simulations
    def genid = Channel.of((1..params.ngens)) 
    def runid = Channel.of((1..params.nruns))


    def nhosts = params.nhosts
    def nsymbionts = params.nsymbionts

    int niter = (int)params.niter
    int freq_subsample = (int)params.subsample

    // Generate data from a coalescent model
    generate_trees_and_interactions(
        genid,
        nsymbionts,
        nhosts,
    )

    // Pass the symbiont tree through revbayes to give it node labels
    rev_annotate_tree(
        genid,
        generate_trees_and_interactions.out.symbiont_tree.map {gid, tree -> tree}
    )
    // Generate the phyJSON
    generate_phyjson(
        genid,
        rev_annotate_tree.out.rev_tree.map  {gid, tree -> tree},
        generate_trees_and_interactions.out.host_tree.map {gid, tree -> tree},
        generate_trees_and_interactions.out.interactions_csv.map {gid, file -> file}
    )

    clean_phyjson(
        genid,
        generate_phyjson.out.dirty_phyjson.map {gid, tree -> tree}
    )
    if (params.run_treeppl) {

        def model_flag_combs
        if (!params.models) {
            // Create all the combinations of compiler flags we desire
            // and give each a unique id
            // If we want to run several models we can specify their keys as --treeppl_model_name <model_key1>,<model_key2> etc.
            def model_names = Channel.fromList(params.treeppl_model_name.tokenize(','))
            def model_dir = Channel.of(params.treeppl_model_dir)
            def inference_flags
            if (params.inference_algorithm == "mcmc-lw-dk") {

                inf_str = "-m mcmc-lw-dk --align --cps none --kernel"
                // Create all combinations of drift we want
                def drift_scale = Channel.of(1.0, 0.1)
                def gprob = Channel.of(0.0)
                inference_flags = drift_scale.combine(gprob)
                    .map {ds, gp -> "${inf_str} --drift ${ds} --mcmc-lw-gprob ${gp}"}

            } else if (params.inference_algorithm in ["smc-apf", "smc-bpf", "mcmc-naive"]) {

                // For these algorithms we don't have additional flags to specify
                inf_str = "-m ${params.inference_algorithm}"
                inference_flags = Channel.of(inf_str).first()

            } else {

                // Use SMC-BPF as default
                inf_str = "-m smc-bpf"
                inference_flags = Channel.of(inf_str).first()

            }

            // Create a channel model and flag combinations
            model_flag_combinations = model_dir.combine(model_names).combine(inference_flags)

        } else {
            // We have a model specification file, use that instead
            model_flag_combinations = Channel.fromPath(params.models)
                | splitCsv(sep:"\t", header:["model_dir", "model_name", "inference_flags"])
                | map {row -> [row.model_dir, row.model_name, row.inference_flags]}
        }
        // Add the rest of the parameters, and the compile id
        compile_id = 0
        compile_config_ch = runid.combine(model_flag_combinations)
            .map {rid, md, mn, flags -> [compile_id++, rid, md, mn, flags]}

        // Save the compile configuration corresponding to each compile id to a file
        compile_config_ch.collectFile(
            name: "compile_id_to_configuration.csv",
            storeDir: file(params.bindir),
            newLine: true
        ) {cid, rid, md, mn, flags -> "$cid\t$rid\t$md\t$mn\t${flags}"}

        // Create all binaries we require
        compile_in_ch = compile_config_ch
            .map {
                cid, rid, md, mn, flags ->
                [cid, rid, md, mn, "$baseDir/models/${params.version_model_dir}/${md}/${mn}.tppl", flags]
            }
        compile_model(compile_in_ch)

        // Create the in file channel
        // -- all combinations of compiler flags and data generations
        treeppl_in_ch = compile_model.out.hostrep_bin.combine(
            clean_phyjson.out.phyjson
        )

        // Run the treeppl implementation
        treeppl_out_ch = run_hostrep_treeppl(
            treeppl_in_ch,
            niter,
        ) 
    }

    if (params.run_revbayes) {
        rev_bayes_in_ch = runid.combine(
            generate_trees_and_interactions.out.symbiont_tree
            .join(generate_trees_and_interactions.out.host_tree)
            .join(generate_trees_and_interactions.out.interactions_nex)
        )

        revbayes_out_ch = run_hostrep_revbayes(
            rev_bayes_in_ch,
            niter,
            freq_subsample,
        )
    }
}
