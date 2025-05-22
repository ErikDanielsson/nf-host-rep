from pathlib import Path
import os
import argparse
import pandas as pd

preamble_str = """---
title: "TreePPL host repertoire"
format:
  html:
    code-fold: true
  pdf:
    code-fold: false
    code-tools: false
    echo: false
    fig-width: 6
    fig-height: 4
execute:
  freeze: auto
jupyter: phylogenetics
---
"""

python_imports_cell_template = """
# Add the python helper directory to the path

import sys
import os
import shutil
from pathlib import Path

CWD = Path(os.getcwd())
sys.path.append(str(CWD / "___python_helpers_dir___"))

# Import the python helpers
import proc_output

import arviz as az
import xarray as xr
import matplotlib.pyplot as plt
import warnings
import pandas as pd
import numpy as np

# Suppress warning from arviz
# -- it complains that all samples are the same
warnings.filterwarnings(action="ignore", module=r"arviz")

# Set global script options, the rest will follow from these
RUN_NAME = "___run_name___"
CACHE = True  # Create a cache of the files to speed up next render?
CLEAR_CACHE = False  # Clear the cache from a previous render?

# Check if we should clear the cache
if not CACHE or CLEAR_CACHE:
    proc_output.clear_temp_dir(tempdir_suffix=RUN_NAME)

# Set the names of the output directories
glob_outdir = CWD / "simulations"
outdir = glob_outdir / RUN_NAME
datadir = outdir / "data"
simdir = outdir / "sims"
bindir = outdir / "bins"
compile_params_path = bindir / "compile_id_to_configuration.csv"
data_params_path = datadir / "param_id_to_configuration.csv"
"""

read_outputs_cell_template = """
# Read files
df = proc_output.get_files_in_dir(
    simdir,
    {
        "tppl": proc_output.get_tppl_output_pattern(),
        "rb": proc_output.get_rb_output_pattern(),
    },
)
df_fn = df.copy()

df = proc_output.create_inference_data_df(
    df,
    {
        "tppl": lambda fn: proc_output.___tppl_file_reader___(
            fn, with_file=CACHE, tempdir_suffix=RUN_NAME
        ),
        "rb": lambda fn: proc_output.read_rb_file(
            fn, with_file=CACHE, tempdir_suffix=RUN_NAME
        ),
    },
    ___burnin___,
    ___subsample_period___,
)
dfs = {k: v for k, v in df.groupby("file_type")}

has_rb = "rb" in dfs
has_tppl = "tppl" in dfs
if has_rb:
    df_rb = dfs["rb"]

if has_tppl:
    df_tppl = dfs["tppl"]
"""

compile_params_text = """
### TreePPL compile params
The following TreePPL models were compiled
"""

compile_params_cell_template = """
compile_params = proc_output.parse_compile_params(compile_params_path)
compile_params
"""

data_params_text = """
The models were run on dataset generated with the following parameters
"""

data_params_cell_template = """
data_params = proc_output.parse_data_params(data_params_path)
data_params.sort_values(by="param_id")
"""

missing_simulations_text = """
### Missing simulations
The following simulations failed to finish (too great RAM requirement)
"""
missing_simulations_cell_template = """
missing_df = proc_output.get_missing_params(df, compile_params_path)
if not missing_df.empty:
    print(missing_df)
else:
    print("All runs finished!")
"""

tppl_data_tree_plots_cell_template = """
df_sim = proc_output.get_files_in_dir(
    datadir,
    {
        "simtrees": proc_output.get_simtree_output_pattern(),
    },
    fields=["param_id", "genid"],
)
df_trees = proc_output.get_trees(df_sim, file_type="simtrees")
tree_plot_paths = {}
for _, row in df_trees.iterrows():
    tree_plot_path = proc_output.create_graphviz_tree_plot(
        row.trees, f"{row.genid}.{row.param_id}", RUN_NAME
    )
    tree_plot_paths[(row.genid, row.param_id)] = tree_plot_path
"""

tppl_log_text = """
### Debug output
Below is a summary of the debug available from the MCMC
"""
tppl_log_cell_template = """
log_df = proc_output.get_files_in_dir(
    simdir,
    {
        "tppl_log": proc_output.get_tppl_log_pattern(),
    },
)
parsed_log_df = proc_output.parse_tppl_logs(log_df, compile_params)
columns = [
    "model_dir",
    "model_name",
    "count",
    "avg_accept",
    "avg_duration_ms",
    "param_id",
    "genid",
    "runid",
]
parsed_log_df.sort_values(by=["compile_id", "param_id"])[columns]
"""


create_multi_chain_tppl_cell_template = """
df_tppl_with_compile_params = proc_output.add_compile_params(
    df_tppl, compile_params
)
reduced_df_tppl = proc_output.create_multi_chain_dataset_df(
    df_tppl_with_compile_params, ___groupby_keys___
)
"""
create_multi_chain_rb_cell_template = """
reduced_df_rb = proc_output.create_multi_chain_dataset_df(
    df_rb,
    ___groupby_keys___,
)
"""

rb_trace_plot_text = """
### RevBayes trace plot___title___
This is a reference run with the model implemented in RevBayes
"""

rb_trace_plot_cell_template = """
# | label: fig-___fig_counter___
# | fig-cap: "___nsamples___ samples from the RevBayes implementation"
# | out.width: '0.5\\linewidth'
# | out.height: '80%'
# | fig.align: center
az.plot_trace(
    reduced_df_rb.loc[___groupby_value___, "multi_channel"],
    figsize=(6, 8),
    compact=True
)
plt.tight_layout()
plt.show()
"""

tppl_trace_plot_text = """
### TreePPL trace plot___title___
"""

tppl_trace_plot_cell_template = """
# | label: fig-trace-tppl-___fig_counter___
# | fig-cap: "TreePPL implementation___title___"
# | out.width: '80%'
# | out.height: '80%'
# | fig.align: center
az.plot_trace(
    reduced_df_tppl.loc[___groupby_value___, "multi_channel"],
    figsize=(6, 8),
    compact=True
)
plt.tight_layout()
plt.show()
"""
rb_ESS_text = """
### RevBayes ESS plots
"""
ess_plot_cell_template = """
fig, ax = plt.subplots(figsize=(8, 10))
proc_output.ess_bar_plot(group_df_ess.get_group(___groupby_value___), ax=ax)
fig.suptitle(___title___)
fig.tight_layout()
plt.show()
"""

rhat_cell_template = """
reduced_df = pd.concat([reduced_df_rb.reset_index(), reduced_df_tppl.reset_index()])
rhat_df = proc_output.get_rhat_dataframe(reduced_df)
rhat_df.sort_values(by=["genid", "param_id", "file_type", "model_dir", "model_name"])
"""

tppl_ESS_text = """
### TreePPL ESS plots
"""

ess_cell_template = """
groupby_tppl = ___groupby_keys_tppl___ 
groupby_rb = ___groupby_keys_rb___ 
groupby_data = ___groupby_keys_data___
df_ess = proc_output.get_all_ess_df(df_tppl_with_compile_params, df_rb, groupby_tppl, groupby_rb, groupby_data, tempdir_suffix=RUN_NAME)
group_df_ess = df_ess.groupby(groupby_data)
"""

tppl_tree_plots_cell_template = """
df_trees = proc_output.get_trees(df_fn)
tree_plot_paths = []
for _, row in df_trees.iterrows():
    tree_plot_path = proc_output.create_graphviz_tree_plot(
        row["trees"], f"{row.genid}.{row.param_id}.{row.compile_id}", RUN_NAME
    )
    tree_plot_paths.append(tree_plot_path)
"""


"""
likelihood = (
    reduced_df_tppl.loc[("likelihood", "uniformization_JC_weight", 1), "multi_channel"]
).posterior.likelihood.to_numpy()
print(likelihood)
fig, axs = plt.subplots(2, 2)
axs = axs.flatten()
sns.histplot(likelihood[0, :], ax=axs[0])
sns.histplot(likelihood[1, :], ax=axs[1])
axs[2].plot(range(10001), sorted(likelihood[0, :]))
axs[3].plot(range(10001), sorted(likelihood[1, :]))
print(likelihood[0, :].mean())
print(likelihood[1, :].mean())
plt.show()
"""

FIG_COUNTER = 0


def get_fig_counter():
    global FIG_COUNTER
    FIG_COUNTER += 1
    return str(FIG_COUNTER - 1)


def subst_variables(cell_str, variables={}):
    for indicator, value in variables.items():
        cell_str = cell_str.replace(indicator, str(value))

    print("#### EXECUTING ####")
    print()
    print(cell_str)
    print()
    print("###################")
    print()
    return cell_str


def create_quarto_cell(cell_str, variables={}):
    cell_start = "```{python}\n"
    cell_end = "\n```\n"
    for indicator, value in variables.items():
        cell_str = cell_str.replace(indicator, str(value))
    return cell_start + cell_str + cell_end


def create_text(text_str, variables={}):
    for indicator, value in variables.items():
        text_str = text_str.replace(indicator, str(value))
    return text_str


def create_single_fig(fig_name, image_fn, caption, title=None):
    lines = []
    if title:
        lines.append(title + "\n")
    lines.append(f"::: {{}}\n\n")

    lines.append(
        f"![{caption}]({Path(image_fn).relative_to(Path(os.getcwd()))}){{#fig-{fig_name}-plots fig-pos=H width=80%}}\n\n"
    )
    lines.append(":::\n")
    return "".join(lines)


def create_multi_fig(fig_name, images_and_captions):
    start = f"::: {{#fig-{fig_name}-plots layout-ncol=2}}\n\n"
    lines = [start]
    for image_fn, caption in images_and_captions:
        lines.append(
            f"![{caption}]({Path(image_fn).relative_to(Path(os.getcwd()))})\n\n"
        )
    lines.append(":::\n")
    return "".join(lines)


def trace_title(gb_keys, gb_vals):
    title = ""
    if "model_dir" in gb_keys:
        title += f": mcmc type `{gb_vals[gb_keys.index('model_dir')]}`,"
    if "model_name" in gb_keys:
        title += f" model type `{gb_vals[gb_keys.index('model_name')]}`"
    return title


def generate_report(
    run_name,
    burnin=0,
    subsample_period=1,
    make_tppl_trace=True,
    make_rb_trace=True,
    make_tppl_tree_plot=False,
    debug_output=False,
    incremental=True,
):
    # Set .gen.qmd as file ending so that we can remove
    # any automatically generated report files
    ph_path_from_base_dir = "python_helpers/"
    ph_path_from_self = "."
    report_name = f"report_{run_name}.gen.qmd"

    data_groupby_keys = ["genid", "param_id"]
    groupby_keys_tppl = ["file_type", "model_dir", "model_name"] + data_groupby_keys
    groupby_keys_rb = ["file_type"] + data_groupby_keys

    groupby_tppl_key = lambda x: (x[2], x[1], x[0])
    groupby_rb_key = lambda x: (x[2], x[1])

    # Execute the python code we need to determine what pipeline was executed
    global_variables = {}
    exec(
        subst_variables(
            python_imports_cell_template,
            {
                "___run_name___": run_name,
                "___python_helpers_dir___": ph_path_from_self,
            },
        ),
        global_variables,
    )
    exec(
        subst_variables(
            read_outputs_cell_template,
            {
                "___burnin___": burnin,
                "___subsample_period___": subsample_period,
                "___tppl_file_reader___": (
                    "read_tppl_incremental_file" if incremental else "read_tppl_file"
                ),
            },
        ),
        global_variables,
    )
    exec(
        subst_variables(
            data_params_cell_template,
            {},
        ),
        global_variables,
    )
    if global_variables["has_tppl"]:
        exec(
            subst_variables(compile_params_cell_template),
            global_variables,
        )
        exec(
            subst_variables(
                create_multi_chain_tppl_cell_template,
                {"___groupby_keys___": groupby_keys_tppl},
            ),
            global_variables,
        )
        if debug_output:
            exec(
                subst_variables(
                    tppl_log_cell_template,
                    {},
                ),
                global_variables,
            )
        if make_tppl_tree_plot:
            exec(
                subst_variables(tppl_tree_plots_cell_template, {}),
                global_variables,
            )
        df_groupby_values_tppl = global_variables["reduced_df_tppl"].index.to_frame(
            index=False
        )
        df_groupby_values = df_groupby_values_tppl
        # groupby_values_tppl.sort(key=groupby_tppl_key)
    if global_variables["has_rb"]:
        exec(
            subst_variables(
                create_multi_chain_rb_cell_template,
                {"___groupby_keys___": groupby_keys_rb},
            ),
            global_variables,
        )
        df_groupby_values_rb = global_variables["reduced_df_rb"].index.to_frame(
            index=False
        )
        if global_variables["has_tppl"]:
            df_groupby_values = pd.concat(
                [df_groupby_values_tppl, df_groupby_values_rb], ignore_index=True
            )
        else:
            df_groupby_values = df_groupby_values_rb
    exec(
        subst_variables(
            tppl_data_tree_plots_cell_template,
            {},
        ),
        global_variables,
    )
    if global_variables["has_rb"] and global_variables["has_tppl"]:
        exec(
            subst_variables(
                ess_cell_template,
                {
                    "___groupby_keys_tppl___": groupby_keys_tppl,
                    "___groupby_keys_rb___": groupby_keys_rb,
                    "___groupby_keys_data___": data_groupby_keys,
                },
            ),
            global_variables,
        )

    # Generate the report
    with open(report_name, "w") as fh:
        fh.write(preamble_str)
        fh.write(
            create_quarto_cell(
                python_imports_cell_template,
                {
                    "___run_name___": run_name,
                    "___python_helpers_dir___": ph_path_from_base_dir,
                },
            )
        )
        fh.write(
            create_quarto_cell(
                read_outputs_cell_template,
                {
                    "___burnin___": burnin,
                    "___subsample_period___": subsample_period,
                    "___tppl_file_reader___": (
                        "read_tppl_incremental_file"
                        if incremental
                        else "read_tppl_file"
                    ),
                },
            )
        )
        fh.write(create_text(data_params_text))
        fh.write(create_quarto_cell(data_params_cell_template))

        if global_variables["has_tppl"]:
            fh.write(create_text(compile_params_text))
            fh.write(create_quarto_cell(compile_params_cell_template))
            fh.write(create_text(missing_simulations_text))
            fh.write(create_quarto_cell(missing_simulations_cell_template))
            if debug_output:
                fh.write(tppl_log_text)
                fh.write(create_quarto_cell(tppl_log_cell_template))

        if global_variables["has_rb"]:
            fh.write(
                create_quarto_cell(
                    create_multi_chain_rb_cell_template,
                    {
                        "___groupby_keys___": groupby_keys_rb,
                    },
                )
            )
        if global_variables["has_tppl"]:
            fh.write(
                create_quarto_cell(
                    create_multi_chain_tppl_cell_template,
                    {"___groupby_keys___": groupby_keys_tppl},
                )
            )

        if global_variables["has_rb"] and global_variables["has_tppl"]:
            fh.write(
                create_quarto_cell(
                    ess_cell_template,
                    {
                        "___groupby_keys_tppl___": groupby_keys_tppl,
                        "___groupby_keys_rb___": groupby_keys_rb,
                        "___groupby_keys_data___": data_groupby_keys,
                    },
                )
            )
            fh.write(
                create_quarto_cell(
                    rhat_cell_template,
                    {},
                )
            )

        fh.write(
            create_quarto_cell(
                tppl_data_tree_plots_cell_template,
                {},
            )
        )
        data_tree_plot_paths = global_variables["tree_plot_paths"]
        data_params = global_variables["data_params"]
        data_params_mi = data_params.set_index("param_id")
        plots_and_titles = {}
        for genid, param_id in sorted(data_tree_plot_paths.keys()):
            param_row = data_params_mi.loc[param_id, :]
            title = (
                f"Tree {genid}, history with params: mu = {param_row.mu},"
                f" beta = {param_row.beta},"
                f" lambda01 = {param_row.lambda01},"
                f" lambda10 = {param_row.lambda10},"
                f" lambda12 = {param_row.lambda12},"
                f" lambda21 = {param_row.lambda21},"
            )
            path = data_tree_plot_paths[(genid, param_id)]
            plots_and_titles[(genid, param_id)] = (path, title)
        if not (global_variables["has_tppl"] and make_tppl_trace) and not (
            global_variables["has_rb"] and make_rb_trace
        ):
            fh.write(create_multi_fig("data_tree_plots", plots_and_titles.values()))

        if True:
            df_groupby_values = df_groupby_values.sort_values(
                by=["genid", "param_id", "file_type"]  # , "model_dir", "model_name"]
            )
            grouped = df_groupby_values.groupby(by=["genid", "param_id"])

            for name, group in grouped:
                tree_plot_fn, tree_plot_caption = plots_and_titles[name]
                fh.write(
                    create_single_fig(
                        f"tree-fig-{int(name[0])}-{int(name[1])}",
                        image_fn=tree_plot_fn,
                        caption=tree_plot_caption,
                        title=f"# Run with tree {int(name[0])} with param set {int(name[1])}",
                    )
                )
                fh.write(
                    create_quarto_cell(
                        ess_plot_cell_template,
                        {
                            "___groupby_value___": name,
                            "___title___": f"'ESS for runs with tree {int(name[0])} and param set {int(name[1])}'",
                        },
                    )
                )

                for _, row in group.iterrows():
                    if row.file_type == "rb":
                        groupby_value = (
                            row.file_type,
                            row.genid,
                            row.param_id,
                        )
                        """
                        fh.write(
                            create_text(
                                rb_trace_plot_text,
                                {
                                    "___title___": trace_title(
                                        groupby_keys_rb, groupby_value
                                    )
                                },
                            )
                        )
                        """
                        fh.write(
                            create_quarto_cell(
                                rb_trace_plot_cell_template,
                                {
                                    "___fig_counter___": get_fig_counter(),
                                    "___nsamples___": 2500,
                                    "___groupby_value___": groupby_value,
                                },
                            )
                        )
                    else:
                        groupby_value = (
                            row.file_type,
                            row.model_dir,
                            row.model_name,
                            row.genid,
                            row.param_id,
                        )
                        """
                        fh.write(
                            create_text(
                                tppl_trace_plot_text,
                                {
                                    "___title___": trace_title(
                                        groupby_keys_tppl, groupby_value
                                    )
                                },
                            )
                        )
                        """
                        fh.write(
                            create_quarto_cell(
                                tppl_trace_plot_cell_template,
                                {
                                    "___fig_counter___": get_fig_counter(),
                                    "___nsamples___": str(2500),
                                    "___groupby_value___": groupby_value,
                                    "___groupby_keys___": groupby_keys_tppl,
                                    "___title___": trace_title(
                                        groupby_keys_tppl, groupby_value
                                    ),
                                },
                            )
                        )

        if global_variables["has_tppl"] and make_tppl_tree_plot:
            fh.write(
                create_quarto_cell(
                    tppl_tree_plots_cell_template,
                    {},
                )
            )
            tree_plot_paths = global_variables["tree_plot_paths"]
            fh.write(
                create_multi_fig(
                    "tree_plots",
                    [(fn, str(k)) for k, fn in tree_plot_paths.items()],
                )
            )
        """
        if global_variables["has_rb"]:
            fh.write(
                create_text(
                    rb_ESS_text,
                    {},
                )
            )
            fh.write(
                create_quarto_cell(
                    rb_ESS_plot_cell_template,
                    {},
                )
            )
        if global_variables["has_tppl"]:
            fh.write(
                create_text(
                    tppl_ESS_text,
                    {},
                )
            )
            fh.write(
                create_quarto_cell(
                    tppl_ESS_plots_cell_template,
                    {"___groupby_keys___": groupby_keys_tppl},
                )
            )
        """


def main():

    parser = argparse.ArgumentParser(
        prog="nf-host-rep-report-generator",
        description="Generate a report for a run of the nf-host-rep pipeline",
        epilog="",
    )
    parser.add_argument("run-name")
    parser.add_argument(
        "--no-tppl-trace",
        action="store_true",
        help="Generate trace plots for TreePPL output",
    )
    parser.add_argument(
        "--no-rb-trace",
        action="store_true",
        help="Generate trace plots for RevBayes output",
    )
    parser.add_argument(
        "--tppl-tree",
        action="store_true",
        help="Generate tree plots for TreePPL output",
    )
    parser.add_argument(
        "--debug-output",
        type=bool,
        default=False,
        help="Generate a summary of the debug output from TreePPL",
    )
    parser.add_argument(
        "--incremental-printing",
        type=bool,
        default=True,
        help="Indicate whether incremental printing was used by TreePPL",
    )
    parser.add_argument(
        "--burnin",
        type=int,
        default=0,
        help="Indicate what burnin period to use",
    )
    parser.add_argument(
        "--subsample-period",
        type=int,
        default=1,
        help="Indicate what subsampling to use (on top of subsampling from infernence methods)",
    )
    args = parser.parse_args()
    run_name = getattr(args, "run-name")
    generate_report(
        run_name,
        make_tppl_trace=not args.no_tppl_trace,
        make_rb_trace=not args.no_rb_trace,
        make_tppl_tree_plot=args.tppl_tree,
        incremental=args.incremental_printing,
        debug_output=args.debug_output,
        burnin=args.burnin,
        subsample_period=args.subsample_period,
    )


if __name__ == "__main__":
    main()
