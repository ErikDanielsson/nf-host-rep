nextflow.enable.dsl=2

include {
    generate_trees;
    clean_phyjson;
    tree_phyjson;
} from "./modules/gendata/trees"

include {
    rev_annotate_tree as rev_annotate_tree_host;
    rev_annotate_tree as rev_annotate_tree_symbiont;
} from "./modules/gendata/trees"

include {
    generate_random_interactions;
} from "./modules/gendata/interactions/random"

include {
    compile_model;
    run_hostrep_treeppl;
} from "./modules/treeppl"

include {
    run_hostrep_revbayes
} from "./modules/revbayes"


include {
    revbayes_interactions;
    clean_rb_csv;
} from "./modules/gendata/interactions/revbayes"

include {
    compile_interactions_tppl;
    run_interactions_tppl;
    add_params_phyjson;
    interactions_json_to_csv;
    interactions_csv_to_nex;
} from "./modules/gendata/interactions/treeppl"

include {
    add_interactions_to_phyjson
} from "./modules/gendata/interactions/helpers"

workflow {
    // Define the simulations
    def runid = Channel.of((1..params.nruns))

    def nhosts = params.nhosts
    def nsymbionts = params.nsymbionts

    int niter = (int)params.niter

    tppl_lib_ch = Channel.fromPath("${params.tppl_lib_path}/*").collect()
    param_id = 0
    params_config_ch = Channel.fromPath("$baseDir/${params.params_config}")
        | splitCsv(sep:"\t", header: true)
        | map {row -> [
            param_id++,
            row.mu,
            row.beta,
            row.lambda01,
            row.lambda10,
            row.lambda12,
            row.lambda21,
        ]}

    params_config_ch.collectFile(
        name: "param_id_to_configuration.csv",
        storeDir: file(params.datadir),
        newLine: true
    ) {pid, mu, beta, l0, l1, l2, l3 -> "$pid\t$mu\t$beta\t$l0\t$l1\t$l2\t$l3"}



    def host_tree_ch
    def symbiont_tree_ch
    def genid
    if (params.empirical_trees) {
        def tree_id = 0
        def empirical_trees_ch = Channel.fromPath("${params.empirical_trees}/*", type: 'dir')
            .map { dir -> 
                def h = dir.listFiles().find { it.name ==~  /.*\.host\..*\.phy/}
                def s = dir.listFiles().find { it.name ==~ /.*\.symbiont\..*\.phy/}
                return (h && s) ? tuple(h, s) : null
            }
            .filter {it != null}
            .map {s, h -> [tree_id++, s, h]}
        host_tree_ch = empirical_trees_ch.map { tid, hpath, spath -> [tid, hpath]}
        symbiont_tree_ch = empirical_trees_ch.map { tid, hpath, spath -> [tid, spath]}
        genid = empirical_trees_ch.map { tid, hpath, spath -> tid }
    } else {
        genid = Channel.of((1..params.ngens)) 
        // Generate data from a coalescent model
        generate_trees(
            genid,
            nsymbionts,
            nhosts
        )
        symbiont_tree_ch = generate_trees.out.symbiont_tree
        host_tree_ch = generate_trees.out.host_tree
    }

    // Pass the symbiont tree through revbayes to give it node labels
    rev_annotate_tree_symbiont(
        symbiont_tree_ch
    )

    rev_annotate_tree_host(
        host_tree_ch
    )

    // Construct a phyjson file containing the symbiont tree
    // and the distance metric induced by the host tree
    // The cleaning step is necessary due to the fact that 
    // R does not distinguish floats and ints
    partial_phyjson_ch = tree_phyjson(
        rev_annotate_tree_symbiont.out.rev_tree.join(
            host_tree_ch
        )
    ) | clean_phyjson

    def interactions_csv_ch
    def interactions_nex_ch
    if (params.interactions == "revbayes") {

        revbayes_interactions_in_ch = rev_annotate_tree_symbiont.out.rev_tree.join(
           host_tree_ch 
        ).join(
            params_config_ch
        )

        revbayes_interactions(
            revbayes_interactions_in_ch
        )

        clean_rb_csv(revbayes_interactions.out.interactions_csv)

        interactions_nex_ch = revbayes_interactions.out.interactions_nex
        interactions_csv_ch = clean_rb_csv.out.interactions_csv
        
    } else if (params.interactions == "treeppl") {
        tppl_sim_ch = params_config_ch.map { row -> row[0] }.combine(Channel.of("$baseDir/bin/simulate.tppl").first())
        compile_interactions_tppl(
            tppl_sim_ch,
            "-m mcmc-lightweight --align --cps full --kernel --sampling-period 1 --incremental-printing --debug-iterations",
        )
        pid = 0
        params_in_ch = partial_phyjson_ch
            .combine(params_config_ch)
            .map {
                gid, pjs, pid, m, b, l1, l2, l3, l4 ->
                [pid, gid, pjs, m, b, l1, l2, l3, l4]
            }
        add_params_phyjson(params_in_ch)

        // Join the channels by the parameter id
        interactions_sim_in_ch = compile_interactions_tppl.out.sim_bin.combine(
            add_params_phyjson.out.phyjson, by: 0
        )

        run_interactions_tppl(
            interactions_sim_in_ch
        )

        interactions_csv_ch = interactions_json_to_csv(
            run_interactions_tppl.out.interactions_json
            .combine(
                rev_annotate_tree_host.out.name_map, by: 0
            ).combine(
                rev_annotate_tree_symbiont.out.name_map, by: 0
            )
        )

        interactions_nex_ch = interactions_csv_to_nex(
            interactions_csv_ch
        )

    } else {
        /* 
         * This will generate a dataset without regards to the phylogenetic signal
         */ 
        generate_random_interactions(
            genid,
            nsymbionts,
            nhosts
        )
        interactions_csv_ch = generate_random_interactions.out.interactions_csv
        interactions_nex_ch = generate_random_interactions.out.interactions_nex

    }

    add_interactions_to_phyjson(
        partial_phyjson_ch.join(
            rev_annotate_tree_symbiont.out.name_map
        ).combine(
            interactions_csv_ch, by: 0
        )
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
        compile_model(compile_in_ch, tppl_lib_ch, niter)

        // Create the in file channel
        // -- all combinations of compiler flags and data generations
        treeppl_in_ch = compile_model.out.hostrep_bin.combine(
            add_interactions_to_phyjson.out.phyjson
        )

        // Run the treeppl implementation
        treeppl_out_ch = run_hostrep_treeppl(
            treeppl_in_ch,
            niter,
        ) 
    }

    if (params.run_revbayes) {
        rev_bayes_in_ch = runid.combine(
           symbiont_tree_ch 
            .join(host_tree_ch)
            .combine(interactions_nex_ch, by: 0)
        )

        revbayes_out_ch = run_hostrep_revbayes(
            rev_bayes_in_ch,
            niter
        )
    }
}
