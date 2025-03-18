import pandas as pd
import numpy as np
import arviz as az
import xarray as xr
import matplotlib.pyplot as plt
import json
import re
import os
import shlex
import shutil
import logging
from functools import reduce
from pathlib import Path

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


def get_temp_dir(tempdir_suffix=""):
    base_tmp_dir = Path(os.getcwd()) / "quarto_tempdirs"
    if not base_tmp_dir.exists():
        base_tmp_dir.mkdir()
    tempdir_prefix = "qmd_temp"
    tempdir_dirname = tempdir_prefix + "_" + tempdir_suffix
    temp_dir = base_tmp_dir / tempdir_dirname
    if not temp_dir.exists():
        temp_dir.mkdir()
    return temp_dir


def clear_temp_dir(tempdir_prefix="qmd_temp", tempdir_suffix=""):
    tempdir_dirname = tempdir_prefix + "_" + tempdir_suffix
    tempdir_path = Path(os.getcwd()) / tempdir_dirname
    if tempdir_path.exists():
        logger.info(f"Clearing cache at {tempdir_path}")
        shutil.rmtree(tempdir_path)
        return True
    return False


def get_temp_file(fn, tempdir_suffix=""):
    return get_temp_dir(tempdir_suffix=tempdir_suffix) / (Path(fn).stem + ".csv")


def get_tppl_output_pattern():
    return re.compile(r"output\.(\d+)\.(\d+)\.json")


def get_rb_output_pattern():
    return re.compile(r"out\.(\d+)\.(\d+)\.log$")


def get_files_in_dir(
    dir: Path,
    patterns: dict[str, re.Pattern],
    header=["file_type", "genid", "compile_id", "filename"],
):
    fns = {k: {} for k in patterns}
    ll = []
    for file in dir.glob("**/*"):
        fn = file.name
        for k, pattern in patterns.items():
            m = pattern.match(fn)
            if m is not None:
                file_ids = map(int, m.groups())
                genid, compile_id = m.groups(1)
                genid, compile_id = int(genid), int(compile_id)
                if genid not in fns[k]:
                    fns[k][genid] = {}
                ll.append((k, *file_ids, file))
                fns[k][genid][compile_id] = file
    df = pd.DataFrame(ll, columns=header)
    return df, fns


def parse_compile_params(
    fn,
    header=["compile_id", "runid", "model_dir", "model_name", "flags"],
    dtypes=[int, int, str, str],
    parse_flags=True,
):
    df = pd.read_csv(
        fn,
        names=header,
        sep="\t",
        index_col=False,
        dtype={h: d for h, d in zip(header, dtypes)},
    )
    if parse_flags:
        flag_df = parse_sim_flags(df)
        return pd.concat([df, flag_df], axis=1).drop(columns=["flags"])

    return df.replace({np.nan: None})


def parse_sim_flags(df, flag_col="flags"):
    return df[flag_col].apply(parse_cmd_line_flag).apply(pd.Series)


def parse_cmd_line_flag(flag):
    flag_pattern = re.compile(r"--?([\w-]+)")  # Matches the flags
    split_flag = shlex.split(flag)
    flag_idx = [flag_pattern.match(arg_part) for arg_part in split_flag]
    args = {}
    i = 0
    while i < len(split_flag):
        if flag_idx[i]:
            flag_name = flag_idx[i].group(1)
            i += 1
            if i < len(split_flag) and not flag_idx[i]:
                args[flag_name] = split_flag[i]
                i += 1
            else:
                args[flag_name] = True
        else:
            print(split_flag[i])
            i += 1
    return args
    # The flag should always have zero or one arguments!


def add_compile_params(file_df, compile_params_df):
    return file_df.merge(compile_params_df, on="compile_id", how="left")


def get_missing_params(file_df, compile_param_fn):
    compile_params_df = parse_compile_params(compile_param_fn)
    found_compile_id = file_df[["genid", "compile_id"]]
    compile_ids = set(compile_params_df["compile_id"])
    missing_compile_ids = {}
    for genid, group in found_compile_id.groupby("genid"):
        existing_compile_ids = set(group["compile_id"])
        missing_compile_ids[genid] = compile_ids - existing_compile_ids
    missing_df = pd.DataFrame(
        [(gen, cid) for gen, cids in missing_compile_ids.items() for cid in cids],
        columns=["genid", "compile_id"],
    )
    return missing_df.merge(compile_params_df, on="compile_id", how="left")


def read_rb_file(fn, rename=True, with_file=True, tempdir_suffix=""):
    if with_file:
        temp_fn = get_temp_file(fn, tempdir_suffix=tempdir_suffix)
        if temp_fn.exists():
            return pd.read_csv(temp_fn, index_col=0)
    samples = pd.read_csv(fn, sep="\t")
    if rename:
        # Just use the columns we are interested in
        samples = samples[
            [
                "clock_host",
                "phy_scale[1]",
                "switch_rate_0_to_1",
                "switch_rate_1_to_0",
                "switch_rate_1_to_2",
                "switch_rate_2_to_1",
            ]
        ]
        # Rename them
        name_map = {
            "clock_host": "mu",
            "phy_scale[1]": "beta",
            "switch_rate_0_to_1": "lambda_01",
            "switch_rate_1_to_0": "lambda_10",
            "switch_rate_1_to_2": "lambda_12",
            "switch_rate_2_to_1": "lambda_21",
        }
        samples = samples.rename(columns=name_map)
    if with_file:
        samples.to_csv(temp_fn)
    return samples


def read_tppl_file(fn, with_file=True, tempdir_suffix=""):
    if with_file:
        temp_fn = get_temp_file(fn, tempdir_suffix=tempdir_suffix)
        if temp_fn.exists():
            return pd.read_csv(temp_fn, index_col=0)
    with open(fn) as fh:
        parsed_file = json.load(fh)
    rename_lambda = [
        "lambda_01",
        "lambda_10",
        "lambda_12",
        "lambda_21",
    ]

    def extract_params(data_entry):
        # Remove the tree
        data_entry.pop("tree")

        # Flatten the lambda entry
        for i, lambda_val in enumerate(data_entry["lambda"]):
            data_entry[rename_lambda[i]] = lambda_val
        data_entry.pop("lambda")

        return data_entry

    # Process the samples, and remove unnecessary params
    df = pd.DataFrame.from_dict(
        [extract_params(s["__data__"]) for s in parsed_file["samples"]]
    )
    df["weights"] = parsed_file["weights"]
    if with_file:
        df.to_csv(temp_fn)

    return df


def inference_data_from_dataframe(df, chain=0, burnin=0, subsample=1):
    index = df.index[burnin::subsample]
    df = df.iloc[index, :]
    df.loc[:, "chain"] = chain
    df.loc[:, "draw"] = index
    df = df.set_index(["chain", "draw"])
    xdata = xr.Dataset.from_dataframe(df)
    return az.InferenceData(posterior=xdata)


def np_data_from_dataframe(df, chain=0, burnin=0, subsample=1):
    index = df.index[burnin::subsample]
    df = df.iloc[index, :]
    df.loc[:, "chain"] = chain
    df.loc[:, "draw"] = index
    df = df.set_index(["chain", "draw"])
    return df.to_numpy
    return az.InferenceData(posterior=xdata)


def create_inference_data_df(file_df, read_funcs, burnin, subsample=1):
    file_df["inference_data"] = [
        inference_data_from_dataframe(
            read_funcs[row["file_type"]](row["filename"]),
            chain=hash((row["genid"], row["compile_id"])),
            burnin=burnin,
            subsample=subsample,
        )
        for _, row in file_df.iterrows()
    ]
    return file_df


def create_np_data_df(file_df, read_funcs, burnin, subsample=1):
    file_df["inference_data"] = [
        np_data_from_dataframe(
            read_funcs[row["file_type"]](row["filename"]),
            chain=hash((row["genid"], row["compile_id"])),
            burnin=burnin,
            subsample=subsample,
        )
        for _, row in file_df.iterrows()
    ]
    return file_df


def create_multi_chain_dataset_df(
    inference_datas, save_cols, data_col="inference_data"
):
    def reduce_col(inference_col):
        concat = lambda x, y: az.concat(x, y, dim="chain")
        return reduce(concat, inference_col)

    return inference_datas.groupby(save_cols).agg(
        multi_channel=pd.NamedAgg(column=data_col, aggfunc=reduce_col)
    )


def create_multi_chain_dataset(inference_datas, type="genid"):
    concat = lambda x, y: az.concat(x, y, dim="chain")
    if type == "all":
        return reduce(
            concat,
            (
                data
                for genid, run_datas in inference_datas.items()
                for runid, data in run_datas.items()
            ),
        )
    if type == "genid":
        return {
            genid: reduce(concat, run_datas.values())
            for genid, run_datas in inference_datas.items()
        }


def calc_ess_all(datas):
    def xarray_to_dict(xarray):
        data_vars = xarray.to_dict()["data_vars"]
        return {k: v["data"] for k, v in data_vars.items()}

    esses = {
        (genid, runid): xarray_to_dict(az.ess(data, method="mean"))
        for genid, run_datas in datas.items()
        for runid, data in run_datas.items()
    }
    df = pd.DataFrame.from_dict(esses, orient="index")
    return df


def ess_bar_plot(
    ess_df: pd.DataFrame, fig=None, axs=None, c="b", label=None, width=0.9, mv=0.0
):
    n_gen = len(ess_df.index.levels[0])
    n_runs = len(ess_df.index.levels[1])
    if fig is None or axs is None:
        fig, axs = plt.subplots(n_gen, 1)
    ess_df_avg = ess_df.groupby(level=0).mean()
    for i, genid in enumerate(ess_df_avg.index):
        ax = axs[i]
        ax.bar(
            np.arange(len(ess_df_avg.columns)) + mv,
            ess_df_avg.loc[genid, :],
            color=c,
            width=width,
            label=label,
        )
        ax.set_ylabel(f"Genid {genid}")
        ax.set_xticks(range(len(ess_df_avg.columns)), ess_df_avg.columns.to_list())
        ax.legend()
        # ax.set_xticks(ess_df_avg.columns)
    fig.suptitle(f"Average ESS over {n_runs} runs")
    return fig, axs


def ess_box_plot(ess_df):
    n_vars = len(ess_df.columns)
    n_gen = len(ess_df.index.levels[0])
    fig, axs = plt.subplots(n_gen, n_vars)
    for i, genid in enumerate(ess_df.index.levels[0]):
        for j, var_name in enumerate(ess_df.columns):
            ax = axs[i, j]
            ax.boxplot(ess_df.xs(genid, level=0)[var_name])
    return fig, axs


def get_time_files(outdir: Path):
    # Define the regexes for the two types of filenames
    tppl_pattern = re.compile(r"time\.treeppl\.(\d+)\.(\d+)\.txt")
    rb_pattern = re.compile(r"time\.revbayes\.(\d+)\.(\d+)\.txt")

    fns = get_files_in_dir(outdir, {"tppl": tppl_pattern, "rb": rb_pattern})

    return fns["tppl"], fns["rb"]


def parse_time_files(time_files):
    return {
        genid: {runid: proc_time_txt(fn) for runid, fn in run_files.items()}
        for genid, run_files in time_files.items()
    }


def proc_time_txt(fn, type="user"):

    def parse_time_str(time_str):
        """
        Extremely hacky parser for the Swedish time
        string format from bash's time command
        """
        pattern = re.compile(r"(\d+)m(\d+),(\d+)s")
        match = pattern.findall(time_str)
        if len(match) == 0:
            raise Exception("Could not parse time str")
        minutes, sec, frac = match[0]
        return int(minutes) * 60 + int(sec) + int(frac) * 10 ** (-len(frac))

    time_df = pd.read_csv(fn, sep="\t", names=["time"], index_col=0)
    time_dict = time_df.to_dict()["time"]
    time_dict = {type: parse_time_str(time_str) for type, time_str in time_dict.items()}
    return time_dict[type]


def interaction_matrix_pattern():
    return re.compile(r"interactions\.(\d+)\.csv$")


def get_interaction_matrices(data_dir):
    pattern = interaction_matrix_pattern()
    fns = {}
    for file in data_dir.iterdir():
        fn = file.name
        m = pattern.match(fn)
        if m is not None:
            (genid,) = m.groups(1)
            fns[genid] = fn
    return fns
