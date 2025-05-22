import sys
import os
import shutil
from pathlib import Path
from functools import reduce

import arviz as az
import xarray as xr
import scipy
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
import pandas as pd
import numpy as np
import seaborn as sns


def aggr_chains_and_times(
    df_MI: pd.DataFrame,
    aggr_dims=[
        "genid",
        "param_id",
        "file_type",
        "model_dir",
        "model_name",
        "mcmc-lw-gprob",
    ],
    data_col="inference_data",
    time_col="parsed_time",
):
    # Reduces an iterators containing InferenceData into a single InferenceData
    def reduce_col(inference_col):
        def concat(x, y):
            t = az.InferenceData(
                posterior=xr.concat([x.posterior, y.posterior], dim="group")
            )
            return t

        return reduce(concat, inference_col)

    df_aggr = df_MI.groupby(aggr_dims).agg(
        pooled=pd.NamedAgg(column=data_col, aggfunc=reduce_col),
        time_sum=pd.NamedAgg(column=time_col, aggfunc="sum"),
        time_mean=pd.NamedAgg(column=time_col, aggfunc="mean"),
        time_std=pd.NamedAgg(column=time_col, aggfunc="std"),
    )
    print(df_aggr)
    return df_aggr


# Functions for computing convergence and ESS per second
def compute_stats(df):
    def calc_ess(inf_data, method):
        chains = az.InferenceData(
            posterior=inf_data.posterior.rename({"group": "chain", "chain": "group"})
        )
        ess = az.ess(chains, method=method)
        return ess.to_dataframe().squeeze()

    def calc_rhat(inf_data):
        chains = az.InferenceData(
            posterior=inf_data.posterior.rename({"group": "chain", "chain": "group"})
        )

        rhat = az.rhat(chains)
        return rhat.to_dataframe().squeeze()

    def count_groups(inf_data):
        n_chains = inf_data.posterior.coords["group"].size
        return n_chains

    df_bulk_ps = df[["pooled", "time_sum"]].apply(
        lambda r: calc_ess(r["pooled"], "bulk") / r["time_sum"], axis=1
    )
    df_tail_ps = df[["pooled", "time_sum"]].apply(
        lambda r: calc_ess(r["pooled"], "tail") / r["time_sum"], axis=1
    )
    df_bulk = df[["pooled", "time_sum"]].apply(
        lambda r: calc_ess(r["pooled"], "bulk"), axis=1
    )
    df_tail = df[["pooled", "time_sum"]].apply(
        lambda r: calc_ess(r["pooled"], "tail"), axis=1
    )
    df_rhat = df["pooled"].apply(calc_rhat)
    n_chains = df["pooled"].apply(count_groups)
    dfs = (df_bulk_ps, df_tail_ps, df_rhat, df_bulk, df_tail)
    for _df in dfs:
        _df["number_of_chains"] = n_chains
    return dfs


def compute_stats_2(df):
    def calc_ess(inf_data, method):
        ess = az.ess(inf_data, method=method)
        return ess.to_dataframe().squeeze()

    df = df.set_index(
        [
            "genid",
            "param_id",
            "file_type",
            "model_dir",
            "model_name",
            "mcmc-lw-gprob",
            "compile_id",
        ]
    )
    df_bulk = df[["inference_data"]].apply(
        lambda r: calc_ess(r["inference_data"], "bulk"), axis=1
    )
    df_tail = df[["inference_data"]].apply(
        lambda r: calc_ess(r["inference_data"], "tail"), axis=1
    )
    return df_bulk, df_tail


def compute_pairwise_stats(file_type_truth="rb", file_type_comp="tppl"):
    # Compute Kolmogorov-Smirnov and Cramér-von Mises statistic between
    # the reference RevBayes implementation
    pass


def show_df(df, vars):
    df["min"] = df[vars].min(axis=1)
    df["max"] = df[vars].max(axis=1)
    df["median"] = df[vars].median(axis=1)
    df["mean"] = df[vars].mean(axis=1)
    df = df.drop(columns=vars)
    return df.reset_index().sort_values(
        by=[
            "genid",
            "param_id",
            "file_type",
            "model_dir",
            "model_name",
            "mcmc-lw-gprob",
        ]
    )


def conv_df(df, thresh=1.01):
    conv = df[["max", "number_of_chains"]].apply(
        lambda row: row["max"] < thresh, axis=1
    )
    print(conv)
    return (
        df[conv]
        .reset_index()
        .sort_values(
            by=[
                "genid",
                "param_id",
                "file_type",
                "model_dir",
                "model_name",
                "mcmc-lw-gprob",
            ]
        )
    )


def show_non_conv_df(show_rhat_df, thresh=1.01):
    non_conv = show_rhat_df[["max", "number_of_chains"]].apply(
        lambda row: row["max"] > thresh, axis=1
    )
    return (
        show_rhat_df[non_conv]
        .reset_index()
        .sort_values(
            by=[
                "genid",
                "param_id",
                "file_type",
                "model_dir",
                "model_name",
                "mcmc-lw-gprob",
            ]
        )
    )
