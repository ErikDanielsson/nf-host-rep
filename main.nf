nextflow.enable.dsl = 2

include {
    generate_trees ;
    clean_phyjson ;
    tree_phyjson
} from "./modules/gendata/trees"

include {
    rev_annotate_tree as rev_annotate_tree_host ;
    rev_annotate_tree as rev_annotate_tree_symbiont
} from "./modules/gendata/trees"

include {
    generate_random_interactions
} from "./modules/gendata/interactions/random"

include {
    compile_model ;
    run_hostrep_treeppl
} from "./modules/treeppl"

include {
    run_hostrep_revbayes
} from "./modules/revbayes"


include {
    revbayes_interactions ;
    clean_rb_csv
} from "./modules/gendata/interactions/revbayes"

include {
    compile_interactions_tppl ;
    run_interactions_tppl ;
    add_params_phyjson ;
    interactions_json_to_csv ;
    interactions_csv_to_nex
} from "./modules/gendata/interactions/treeppl"

include {
    add_interactions_to_phyjson
} from "./modules/gendata/interactions/helpers"

workflow {
    // Define the simulations
    def runid = Channel.of((1..params.nruns))

    def nhosts = params.nhosts
    def nsymbionts = params.nsymbionts

    def int niter = params.niter

    tppl_lib_ch = Channel.fromPath("${params.tppl_lib_path}/*").collect()

    // Read the simulation parameters from the provided file
    param_id = 0
    params_config_ch = Channel.fromPath("${baseDir}/${params.params_config}")
        | splitCsv(sep: "\t", header: true)
        | map { row ->
            [
                ++param_id,
                row.mu,
                row.beta,
                row.lambda01,
                row.lambda10,
                row.lambda12,
                row.lambda21,
                row.seed,
            ]
        }

    // Save the parameters into the output directory to be read by downstream analyses
    params_config_ch.collectFile(
        name: "param_id_to_configuration.csv",
        storeDir: file(params.datadir),
        newLine: true,
        sort: true,
    ) { pid, mu, beta, l0, l1, l2, l3, s -> "${pid}\t${mu}\t${beta}\t${l0}\t${l1}\t${l2}\t${l3}\t${s}" }

    /*
     * Generate the host and symbiont trees
     * 
     * If empirical trees are provided, we will use those.
     * Otherwise, we will generate trees from a coalescent model.
     * 
     */
    def host_tree_ch
    def symbiont_tree_ch
    def genid
    if (params.empirical_trees) {
        def tree_id = 0
        // Search the provideed direcory for host and symbiont trees
        def empirical_trees_ch = Channel.fromPath("${params.empirical_trees}/*", type: 'dir')
            .map { dir ->
                def h = dir.listFiles().find { it.name ==~ /.*\.host\..*\.phy/ }
                def s = dir.listFiles().find { it.name ==~ /.*\.symbiont\..*\.phy/ }
                return h && s ? tuple(h, s) : null
            }
            .filter { it != null }
            .map { s, h -> [++tree_id, s, h] }
        empirical_trees_ch.collectFile(
            name: "empirical_trees.csv",
            storeDir: file(params.datadir),
            newLine: true,
            sort: true,
        ) { tid, s, h -> "${tid}\t${s}\t${h}" }
        host_tree_ch = empirical_trees_ch.map { tid, hpath, spath -> [tid, hpath] }
        symbiont_tree_ch = empirical_trees_ch.map { tid, hpath, spath -> [tid, spath] }
        genid = empirical_trees_ch.map { tid, hpath, spath -> tid }
    }
    else {
        genid = Channel.of((1..params.ngens))
        // Generate host and symbiont trees from a coalescent model
        generate_trees(
            genid,
            nsymbionts,
            nhosts,
        )
        symbiont_tree_ch = generate_trees.out.symbiont_tree
        host_tree_ch = generate_trees.out.host_tree
    }

    // Pass the trees through RevBayes to give it consistent node labels
    rev_annotate_tree_symbiont(
        symbiont_tree_ch
    )
    rev_annotate_tree_host(
        host_tree_ch
    )

    // Construct a phyjson file containing the symbiont tree
    // and the distance metric induced by the host tree
    // The cleaning step is necessary since R does not distinguish
    // floats and ints
    partial_phyjson_ch = tree_phyjson(
        rev_annotate_tree_symbiont.out.rev_tree.join(
            host_tree_ch
        )
    )
        | clean_phyjson

    /*
     * Generate the interactions
     * 
     * If we use the treeppl script to generate interactions with a phylogenetic signal,
     * we will use the provided parameters in the `params_config_ch` channel.
     * 
     * If we do not use the treeppl script, we will generate interactions
     * without regards to the trees.
     * 
     */
    def interactions_csv_ch
    def interactions_nex_ch
    if (params.interactions == "treeppl") {
        // Create a channel with the parameter id and the seed
        // to be used when compiling the TreePPL script
        tppl_sim_ch = params_config_ch
            .map { pid, mu, beta, l0, l1, l2, l3, s ->
                [pid, s]
            }
            .combine(Channel.of("${baseDir}/bin/simulate_with_subroot.tppl").first())

        // Compile the script
        compile_interactions_tppl(
            tppl_sim_ch,
            "-m mcmc-lightweight --align --cps full --kernel --sampling-period 1 --incremental-printing --debug-iterations",
        )

        // Construct the input phyjson:
        // - the phyjson containing the tree and distance metric
        // - the parameters to be used in the simulation
        params_in_ch = partial_phyjson_ch
            .combine(params_config_ch)
            .map { gid, pjs, pid, m, b, l1, l2, l3, l4, s ->
                [pid, gid, pjs, m, b, l1, l2, l3, l4]
            }
        // Add the parameters to the phyjson
        add_params_phyjson(params_in_ch)

        // Join the channels by the parameter id
        interactions_sim_in_ch = compile_interactions_tppl.out.sim_bin.combine(
            add_params_phyjson.out.phyjson,
            by: 0
        )

        // Generate the interactions
        run_interactions_tppl(
            interactions_sim_in_ch
        )

        // To allow RevBayes to use the interactions we generate 
        // a NeXus via a CSV file.
        interactions_csv_ch = interactions_json_to_csv(
            run_interactions_tppl.out.interactions_json.combine(
                rev_annotate_tree_host.out.name_map,
                by: 0
            ).combine(
                rev_annotate_tree_symbiont.out.name_map,
                by: 0
            )
        )
        interactions_nex_ch = interactions_csv_to_nex(
            interactions_csv_ch
        )
    }
    else {
        // Generate a dataset without phylogenetic signal
        generate_random_interactions(
            genid,
            nsymbionts,
            nhosts,
        )
        interactions_csv_ch = generate_random_interactions.out.interactions_csv
        interactions_nex_ch = generate_random_interactions.out.interactions_nex
    }

    // Add the inteactions produced by either method to the phyjson
    add_interactions_to_phyjson(
        partial_phyjson_ch.join(
            rev_annotate_tree_symbiont.out.name_map
        ).combine(
            interactions_csv_ch,
            by: 0
        )
    )

    /*
     * Run the host-repertoire model (HRM)
     *
     * TreePPL: Run the variants of the HRM specified in the `params.models` file
     * RevBayes: Run the HRM as specified in Braga et al. 2020 paper
     * 
     * Each model will be compiled as many times as specfied with the `--nruns` flag.
     */
    if (params.run_treeppl) {
        // We have a model specification file, use that instead
        model_flag_combinations = Channel.fromPath(params.models)
            | splitCsv(sep: "\t", header: ["model_dir", "model_name", "inference_flags"])
            | map { row -> [row.model_dir, row.model_name, row.inference_flags] }
        // Add the rest of the parameters, and the compile id
        compile_id = 0
        compile_config_ch = runid
            .combine(model_flag_combinations)
            .map { rid, md, mn, flags -> [++compile_id, rid, md, mn, flags] }

        // Save the compile configuration corresponding to each compile id to a file
        compile_config_ch.collectFile(
            name: "compile_id_to_configuration.csv",
            storeDir: file(params.bindir),
            newLine: true,
            sort: true,
        ) { cid, rid, md, mn, flags -> "${cid}\t${rid}\t${md}\t${mn}\t${flags}" }

        // Create all binaries we require
        compile_in_ch = compile_config_ch.map { cid, rid, md, mn, flags ->
            [cid, rid, md, mn, "${baseDir}/models/${params.version_model_dir}/${md}/${mn}.tppl", flags]
        }
        compile_model(compile_in_ch, tppl_lib_ch, niter)

        // Create the in file channel -- all combs of compiler flags and data generations
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
            symbiont_tree_ch.join(host_tree_ch).combine(interactions_nex_ch, by: 0)
        )

        revbayes_out_ch = run_hostrep_revbayes(
            rev_bayes_in_ch,
            niter,
        )
    }
}
