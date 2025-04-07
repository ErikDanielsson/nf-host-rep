from pathlib import Path
import os

preamble_str = """---
title: "TreePPL host repertoire"
format:
  html:
    code-fold: true
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
param_comb_path = bindir / "compile_id_to_configuration.csv"
"""

read_outputs_cell_template = """
# Read files
df, tppl_fns = proc_output.get_files_in_dir(
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
    1,
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
The following TreePPL models were run
"""

compile_params_cell_template = """
compile_params = proc_output.parse_compile_params(param_comb_path)
compile_params
"""

missing_simulations_text = """
### Missing simulations
The following simulations failed to finish (too great RAM requirement)
"""
missing_simulations_cell_template = """
missing_df = proc_output.get_missing_params(df, param_comb_path)
if not missing_df.empty:
    print(missing_df)
else:
    print("All runs finished!")
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
### RevBayes trace plot `___groupby_keys___` = `___groupby_value___`
This is a reference run with the model implemented in RevBayes
"""

rb_trace_plot_cell_template = """
# | label: fig-___fig_counter___
# | fig-cap: "___nsamples___ samples from the RevBayes implementation"
az.plot_trace(
    reduced_df_rb.loc[___groupby_value___, "multi_channel"],
    figsize=(10, 10),
    compact=False
)
plt.show()
"""

tppl_trace_plot_text = """
### TreePPL trace plot `___groupby_keys___` = `___groupby_value___`
"""

tppl_trace_plot_cell_template = """
# | label: fig-trace-tppl-___fig_counter___
# | fig-cap: "___nsamples___ samples from the TPPL implementation  `___groupby_value___`"
az.plot_trace(
    reduced_df_tppl.loc[___groupby_value___, "multi_channel"],
    figsize=(10, 10),
    compact=False
)
plt.show()
"""
rb_ESS_text = """
### RevBayes ESS plots
"""
rb_ESS_plot_cell_template = """
df_ess_rb = proc_output.get_ess_df(df_rb, [])
fig, ax = plt.subplots()
proc_output.ess_bar_plot(df_ess_rb, ax)
fig.suptitle("ESS for the RevBayes run")
fig.set_size_inches(8, 10)
fig.tight_layout()
plt.show()
"""

tppl_ESS_text = """
### TreePPL ESS plots
"""

tppl_ESS_plots_cell_template = """
groupby = ___groupby_keys___ 
df_ess = proc_output.get_ess_df(df_tppl_with_compile_params, groupby)
fig, axs = proc_output.ess_group_bar_plot(df_ess, groupby)
fig.set_size_inches(8, 20)
fig.tight_layout()
plt.show()
"""

tppl_tree_plots_cell_template = """
df_trees = proc_output.get_trees(df_fn)
tree_plot_paths = []
for i in df_trees.index:
    tree_plot_path = proc_output.create_graphviz_tree_plot(
        df_trees.loc[i, "trees"], i, RUN_NAME
    )
    tree_plot_paths.append(tree_plot_path)
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


def create_multi_fig(fig_name, images_and_captions):
    start = f"::: {{#fig-{fig_name}-plots layout-ncol=2}}\n\n"
    lines = [start]
    for image_fn, caption in images_and_captions:
        lines.append(
            f"![{caption}]({Path(image_fn).relative_to(Path(os.getcwd()))})\n\n"
        )
    lines.append(":::\n")
    return "".join(lines)


def generate_report(run_name, burnin=0):
    # Set .gen.qmd as file ending so that we can remove
    # any automatically generated report files
    make_tppl_trace = True
    make_rb_trace = True
    make_tppl_tree_plot = False
    incremental = True
    ph_path_from_base_dir = "python_helpers/"
    ph_path_from_self = "."
    report_name = f"report_{run_name}.gen.qmd"

    groupby_keys_tppl = ["model_dir", "model_name", "genid"]
    groupby_keys_rb = ["file_type", "genid"]

    groupby_tppl_key = lambda x: (x[2], x[1], x[0])

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
                "___tppl_file_reader___": (
                    "read_tppl_incremental_file" if incremental else "read_tppl_file"
                ),
            },
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
        if make_tppl_tree_plot:
            exec(
                subst_variables(tppl_tree_plots_cell_template, {}),
                global_variables,
            )
        groupby_values_tppl = list(global_variables["reduced_df_tppl"].index)
        groupby_values_tppl.sort(key=groupby_tppl_key)
    if global_variables["has_rb"]:
        exec(
            subst_variables(
                create_multi_chain_rb_cell_template,
                {"___groupby_keys___": groupby_keys_rb},
            ),
            global_variables,
        )
        groupby_values_rb = global_variables["reduced_df_rb"].index

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
                    "___tppl_file_reader___": (
                        "read_tppl_incremental_file"
                        if incremental
                        else "read_tppl_file"
                    ),
                },
            )
        )
        if global_variables["has_tppl"]:
            fh.write(create_text(compile_params_text))
            fh.write(create_quarto_cell(compile_params_cell_template))
            fh.write(create_text(missing_simulations_text))
            fh.write(create_quarto_cell(missing_simulations_cell_template))

        if global_variables["has_rb"]:
            fh.write(
                create_quarto_cell(
                    create_multi_chain_rb_cell_template,
                    {
                        "___groupby_keys___": groupby_keys_rb,
                    },
                )
            )

        if global_variables["has_rb"] and make_rb_trace:
            for groupby_value in groupby_values_rb:
                fh.write(
                    create_text(
                        rb_trace_plot_text,
                        {
                            "___groupby_keys___": groupby_keys_rb,
                            "___groupby_value___": groupby_value,
                        },
                    )
                )
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

        if global_variables["has_tppl"]:
            fh.write(
                create_quarto_cell(
                    create_multi_chain_tppl_cell_template,
                    {"___groupby_keys___": groupby_keys_tppl},
                )
            )

        if global_variables["has_tppl"] and make_tppl_trace:
            for groupby_value in groupby_values_tppl:
                fh.write(
                    create_text(
                        tppl_trace_plot_text,
                        {
                            "___groupby_value___": groupby_value,
                            "___groupby_keys___": groupby_keys_tppl,
                        },
                    )
                )
                fh.write(
                    create_quarto_cell(
                        tppl_trace_plot_cell_template,
                        {
                            "___fig_counter___": get_fig_counter(),
                            "___nsamples___": str(2500),
                            "___groupby_value___": groupby_value,
                            "___groupby_keys___": groupby_keys_tppl,
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
                    zip(tree_plot_paths, tree_plot_paths),
                )
            )
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


def main():
    import sys

    run_name = sys.argv[1]
    generate_report(run_name)


if __name__ == "__main__":
    main()
