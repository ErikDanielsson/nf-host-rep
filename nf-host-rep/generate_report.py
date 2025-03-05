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
sys.path.append(str(CWD / "../python_helpers/"))

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
outdir = CWD / RUN_NAME
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

df = proc_output.create_inference_data_df(
    df,
    {
        "tppl": lambda fn: proc_output.read_tppl_file(
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
    df_tppl_with_compile_params, ["___groupby_key___"]
)
"""
create_multi_chain_rb_cell_template = """
reduced_df_rb = proc_output.create_multi_chain_dataset_df(
    df_rb,
    ["___groupby_key___"],
)
"""
rb_trace_plot_text = """
### RevBayes trace plot
This is a reference run with the model implemented in RevBayes
"""

rb_trace_plot_cell_template = """
# | label: fig-___fig_counter___
# | fig-cap: "___nsamples___ samples from the RevBayes implementation"
if has_rb:
    az.plot_trace(reduced_df_rb.loc["__groupby_value___", "multi_channel"], compact=False)
    plt.show()
else:
    print("No RevBayes files found")
"""

tppl_trace_plot_text = """
### TreePPL trace plot ___groupby_key___ = ___groupby_value___
"""

tppl_trace_plot_cell_template = """
# | label: fig-trace-tppl-___fig_counter___
# | fig-cap: "___nsamples___ samples from the TPPL implementation with drift param ___groupby_value___"
az.plot_trace(reduced_df_tppl.loc["___groupby_value___", "multi_channel"], compact=False)
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


def generate_report(run_name, burnin=0):
    # Set .gen.qmd as file ending so that we can remove
    # any automatically generated report files
    report_name = f"report_{run_name}.gen.qmd"

    groupby_key_tppl = "drift"
    groupby_key_rb = "file_type"
    groupby_value_rb = "rb"

    # Execute the python code we need to determine what pipeline was executed
    global_variables = {}
    exec(
        subst_variables(python_imports_cell_template, {"___run_name___": run_name}),
        global_variables,
    )
    exec(
        subst_variables(
            read_outputs_cell_template,
            {"___burnin___": burnin},
        ),
        global_variables,
    )
    exec(
        subst_variables(compile_params_cell_template),
        global_variables,
    )
    exec(
        subst_variables(
            create_multi_chain_tppl_cell_template,
            {"___groupby_key___": groupby_key_tppl},
        ),
        global_variables,
    )
    groupby_values_tppl = global_variables["reduced_df_tppl"].index

    # Generate the report
    with open(report_name, "w") as fh:
        fh.write(preamble_str)
        fh.write(
            create_quarto_cell(
                python_imports_cell_template, {"___run_name___": run_name}
            )
        )
        fh.write(
            create_quarto_cell(
                read_outputs_cell_template,
                {"___burnin___": burnin},
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
                    {"___groupby_key___": groupby_key_rb},
                )
            )
            fh.write(create_text(rb_trace_plot_text))
            fh.write(
                create_quarto_cell(
                    rb_trace_plot_cell_template,
                    {
                        "___fig_counter___": get_fig_counter(),
                        "___nsamples___": 2500,
                        "__groupby_value___": groupby_value_rb,
                    },
                )
            )
        if global_variables["has_tppl"]:
            fh.write(
                create_quarto_cell(
                    create_multi_chain_tppl_cell_template,
                    {"___groupby_key___": groupby_key_tppl},
                )
            )

            for groupby_value in groupby_values_tppl:
                fh.write(
                    create_text(
                        tppl_trace_plot_text,
                        {
                            "___groupby_value___": groupby_value,
                            "___groupby_key___": groupby_key_tppl,
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
                            "___groupby_key___": groupby_key_tppl,
                        },
                    )
                )


def main():
    import sys

    run_name = sys.argv[1]
    generate_report(run_name)


if __name__ == "__main__":
    main()
