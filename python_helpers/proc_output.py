import pandas as pd
import numpy as np
import arviz as az
import graphviz
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
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


def read_tppl_incremental_file(fn, with_file=True, tempdir_suffix=""):
    if with_file:
        temp_fn = get_temp_file(fn, tempdir_suffix=tempdir_suffix)
        if temp_fn.exists():
            return pd.read_csv(temp_fn, index_col=0)

    def extract_params(data_entry):
        rename_lambda = [
            "lambda_01",
            "lambda_10",
            "lambda_12",
            "lambda_21",
        ]
        # Remove the tree
        data_entry.pop("tree")

        # Flatten the lambda entry
        for i, lambda_val in enumerate(data_entry["lambda"]):
            data_entry[rename_lambda[i]] = lambda_val
        data_entry.pop("lambda")

        return data_entry

    samples = []
    with open(fn) as fh:
        for line in fh:
            parsed_line = json.loads(line)
            if "__data__" in parsed_line:
                samples.append(extract_params(parsed_line["__data__"]))

    df = pd.DataFrame.from_dict(samples)
    if with_file:
        df.to_csv(temp_fn)
    return df


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


def get_json_tree_tppl(fn):
    with open(fn) as fh:
        parsed_file = json.load(fh)
    tree_samples = [s["__data__"]["tree"] for s in parsed_file["samples"]]
    return tree_samples


def parse_tree(json_tree_samples):
    def initial_tree(first_tree):
        data = first_tree["__data__"]
        if "right" in data:
            return {
                "label": str(data["label"]),
                "age": data["age"],
                "repertoire": [data["repertoire"]],
                "right": initial_tree(data["right"]),
                "left": initial_tree(data["left"]),
            }
        else:
            return {
                "label": str(data["label"]),
                "age": data["age"],
                "repertoire": [data["repertoire"]],
            }

    def append_repertoires(tree_sample, parsed_tree):
        data = tree_sample["__data__"]
        if "right" in parsed_tree:
            left = parsed_tree["left"]
            right = parsed_tree["right"]
            parsed_tree["repertoire"].append(data["repertoire"])
            append_repertoires(data["left"], left)
            append_repertoires(data["right"], right)
        else:
            parsed_tree["repertoire"].append(data["repertoire"])

    tree = initial_tree(json_tree_samples[0])
    for tree_sample in json_tree_samples[1:]:
        append_repertoires(tree_sample, tree)
    return tree


def tree_map(tree, func, inkeys, outkey):
    if "right" in tree:
        return {
            outkey: func(*[tree[k] for k in inkeys]),
            "right": tree_map(tree["right"], func, inkeys, outkey),
            "left": tree_map(tree["left"], func, inkeys, outkey),
            "label": tree["label"],
            "age": tree["age"],
        }
    else:
        return {
            outkey: func(*[tree[k] for k in inkeys]),
            "label": tree["label"],
            "age": tree["age"],
        }


def tree_avg_rep(tree):
    def rep_to_vec(vec, size=3):
        nvec = np.array(vec)
        return np.eye(size)[nvec]

    def transform_rep_array(rep_array):
        arr = np.array([rep_to_vec(r) for r in rep_array]).mean(axis=0)
        return arr

    return tree_map(tree, transform_rep_array, ["repertoire"], "avg_rep")


heat_map_dir = Path(".heatmaps")


def heat_tree(mat_tree, tempdir):

    def mat_to_heatmap_file(mat, node_name):
        filename = f"heatmap.{node_name}.png"
        file_path = tempdir / filename
        if not file_path.exists():
            plt.figure(figsize=(2, 2), dpi=100)  # Fixed size and resolution
            plt.imshow(mat.T, cmap="binary", aspect="auto")
            # plt.axis("off")  # Hide axes
            plt.xlabel("Host")
            plt.ylabel("State")
            plt.savefig(file_path, bbox_inches="tight", pad_inches=0)
            plt.close()
        return os.path.abspath(str(file_path))

    return tree_map(mat_tree, mat_to_heatmap_file, ["avg_rep", "label"], "hmapfn")


def create_dot_tree(tree, dot):
    if "right" in tree:
        dot.node(
            tree["label"],
            label="",
            image=tree["hmapfn"],
            shape="box",
            width="1.2",
            height="1.2",
            imagescale="true",
        )
        dot.edge(tree["label"], str(tree["right"]["label"]))
        dot.edge(tree["label"], str(tree["left"]["label"]))
        create_dot_tree(tree["left"], dot)
        create_dot_tree(tree["right"], dot)
    else:
        dot.node(
            tree["label"],
            label="",
            image=tree["hmapfn"],
            shape="box",
            width="1.2",
            height="1.2",
            imagescale="true",
        )


def create_graphviz_tree_plot(json_tree, id, tempdir_suffix):
    tempdir = get_temp_dir(tempdir_suffix)
    file_path = tempdir / f"treeplot.{id}"
    if not file_path.exists():
        parsed_tree = parse_tree(json_tree)
        mat_tree = tree_avg_rep(parsed_tree)
        heat_map_tempdir = tempdir / f"heatmaps.{id}"
        heat_map_tempdir.mkdir()
        heat_map_tree = heat_tree(mat_tree, heat_map_tempdir)
        dot = graphviz.Digraph()
        create_dot_tree(heat_map_tree, dot)
        print()
        path = dot.render(file_path, format="png", view=True)
        return path
    else:
        return f"{file_path}.png"


def get_trees(df_fn):
    df_fn = df_fn[df_fn["file_type"] == "tppl"]
    df_fn = df_fn.reset_index()
    df_fn["trees"] = [
        get_json_tree_tppl(row["filename"]) for _, row in df_fn.iterrows()
    ]
    return df_fn


def inference_data_from_dataframe(df, chain=0, burnin=0, subsample=1):
    index = df.index[burnin::subsample]
    df = df.iloc[index, :]
    df.loc[:, "chain"] = chain
    df.loc[:, "draw"] = index
    df = df.set_index(["chain", "draw"])
    xdata = xr.Dataset.from_dataframe(df)
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


def create_multi_chain_dataset_df(
    inference_datas, save_cols, data_col="inference_data"
):
    def reduce_col(inference_col):
        concat = lambda x, y: az.concat(x, y, dim="chain")
        return reduce(concat, inference_col)

    return inference_datas.groupby(save_cols).agg(
        multi_channel=pd.NamedAgg(column=data_col, aggfunc=reduce_col)
    )


def approx_eq(f1, f2):
    APPROX_SENSE = 1e-10
    return abs(f1 - f2) < APPROX_SENSE


def get_ess_df(df, groupby, data_col="inference_data"):

    def xarray_to_series(data: az.InferenceData):
        N = len(data.posterior.draw)
        ess_all = az.ess(data, method="bulk")
        vars = ess_all.keys()
        return pd.Series(
            [float(ess_all[v]) if not approx_eq(ess_all[v], N) else 0 for v in vars],
            index=vars,
        )

    index_name = "index"
    df.index.name = index_name
    df_mi = df.reset_index()
    df_mi = df_mi.set_index(groupby + [index_name])
    df_mi = df_mi[data_col].apply(xarray_to_series)
    df_mi = df_mi.reset_index()
    return df_mi


def ess_group_bar_plot(df_ess, groupby):
    grouped_df_ess = df_ess.groupby(groupby)
    fig, axs = plt.subplots(grouped_df_ess.ngroups, 1, squeeze=False)
    for i, (gname, group) in enumerate(grouped_df_ess):
        ax = axs[i, 0]
        group = group.drop(groupby, axis=1)
        ess_bar_plot(group, ax)
        ax.set_title(f"ESS for model '{': '.join(map(str, gname))}'")
    fig.tight_layout()
    return fig, axs


def calc_ess(df, data_col="inference_data"):
    def xarray_to_series(data: xr.DataArray):
        ess_all = az.ess(data, method="bulk")
        vars = ess_all.keys()
        return pd.Series([float(ess_all[v]) for v in vars], index=vars)

    df_ess = df[data_col].apply(xarray_to_series)
    return df_ess


def ess_bar_plot(df_ess, ax=None):
    df_ess_long = df_ess.melt(id_vars="index", var_name="Variable", value_name="ESS")
    return sns.barplot(
        data=df_ess_long,
        y="Variable",
        x="ESS",
        hue="index",
        legend=False,
        ax=ax,
    )


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
