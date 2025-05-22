import pandas as pd
import numpy as np
import json
import os
import graphviz
import matplotlib.pyplot as plt
from itertools import islice


class RepertoireTree:

    def __init__(self, tree=None, type=None):
        if type == "json":
            self.tree = self.parse_json_tree(tree)
        else:
            self.tree = tree

    def parse_json_tree(self, tree):
        data = tree["__data__"]
        if "right" in data:
            return {
                "label": str(data["label"]),
                "age": data["age"],
                "repertoire": [data["repertoire"]],
                "history_length": [len(data["history"])],
                "right": self.parse_json_tree(data["right"]),
                "left": self.parse_json_tree(data["left"]),
            }
        else:
            return {
                "label": str(data["label"]),
                "age": data["age"],
                "repertoire": [data["repertoire"]],
            }

    def append_json(self, new_tree):
        def append_helper(new_tree, current_tree):
            data = new_tree["__data__"]
            if "right" in current_tree:
                left = current_tree["left"]
                right = current_tree["right"]
                current_tree["repertoire"].append(data["repertoire"])
                current_tree["history_length"].append(len(data["history"]))
                append_helper(data["left"], left)
                append_helper(data["right"], right)
            else:
                current_tree["repertoire"].append(data["repertoire"])

        if self.tree is None:
            self.tree = self.parse_json_tree(new_tree)
        else:
            append_helper(new_tree, self.tree)

    def update_json(self, other_tree):
        def update_helper(their_tree, our_tree):
            data = their_tree["__data__"]
            if "right" in our_tree:
                left = our_tree["left"]
                right = our_tree["right"]
                our_tree["repertoire"] += their_tree["repertoire"]
                our_tree["history_length"] += their_tree["repertoire"]
                update_helper(data["left"], left)
                update_helper(data["right"], right)
            else:
                our_tree["repertoire"] += their_tree["repertoire"]

        update_helper(other_tree, self.tree)

    def map(self, func, inkeys, outkey):
        def map_helper(tree, func, inkeys, outkey):
            if "right" in tree:
                val = (
                    func(*[tree[k] for k in inkeys])
                    if all(map(lambda k: k in tree, inkeys))
                    else None
                )
                return dict(
                    list(tree.items())
                    + [
                        (outkey, val),
                        ("right", map_helper(tree["right"], func, inkeys, outkey)),
                        ("left", map_helper(tree["left"], func, inkeys, outkey)),
                    ]
                )
            else:
                val = (
                    func(*[tree[k] for k in inkeys])
                    if all(map(lambda k: k in tree, inkeys))
                    else None
                )
                return dict(
                    list(tree.items())
                    + [
                        (outkey, val),
                    ]
                )

        new_tree = map_helper(self.tree, func, inkeys, outkey)
        return RepertoireTree(new_tree)

    def filter(self, tree, keys):
        def filter_helper(tree, keys):
            if "right" in tree:
                return dict(
                    [(k, v) for k, v in tree.items() if k in keys]
                    + [
                        ("right", filter_helper(tree["right"], keys)),
                        ("left", filter_helper(tree["left"], keys)),
                    ]
                )
            else:
                return dict([(k, v) for k, v in tree.items() if k in keys])

        new_tree = filter_helper(tree, keys)
        return RepertoireTree(new_tree)

    def create_graphviz_tree(
        self,
        node_args=[lambda t: t["label"]],
        node_kwargs={
            "image": (lambda t: t["hmapfn"]),
            "shape": (lambda t: "box"),
            "width": (lambda t: "1.2"),
            "height": (lambda t: "1.2"),
            "imagescale": (lambda t: "true"),
        },
        edge_args=[
            lambda t, d: t["label"],
            lambda t, d: t[d]["label"],
        ],
        edge_kwargs={
            "headlabel": (lambda t, d: "0"),
            "len": (lambda t, d: str(t["age"] - t[d]["age"])),
        },
    ):
        def create_dot_tree(tree, dot):
            if "right" in tree:
                dot.node(
                    *[arg(tree) for arg in node_args],
                    **{k: v(tree) for k, v in node_kwargs.items()},
                )
                dot.edge(
                    *[arg(tree, "right") for arg in edge_args],
                    **{k: v(tree, "right") for k, v in edge_kwargs.items()},
                )
                dot.edge(
                    *[arg(tree, "left") for arg in edge_args],
                    **{k: v(tree, "left") for k, v in edge_kwargs.items()},
                )
                create_dot_tree(tree["left"], dot)
                create_dot_tree(tree["right"], dot)
            else:
                dot.node(
                    *[arg(tree) for arg in node_args],
                    **{k: v(tree) for k, v in node_kwargs.items()},
                )

        dot = graphviz.Digraph()
        create_dot_tree(self.tree, dot)
        return dot


class PosteriorTrees:
    def __init__(
        self,
        df_fn,
        global_tempdir,
        metadata_cols=[
            "genid",
            "param_id",
            "compile_id",
            "file_type",
            "model_dir",
            "model_name",
        ],
    ):
        self.df_fn = df_fn
        self.df_trees = df_fn[metadata_cols]
        self.tppl_idx = self.df_trees.index[self.df_trees["file_type"] == "tppl"]
        self.rb_idx = self.df_trees.index[self.df_trees["file_type"] == "rb"]
        self.active_idx = pd.Index([])
        self.global_tempdir = global_tempdir

    def get_tppl_trees(self, burnin=0, subsample=1, end=None):
        tppl_idx = self.df_fn["file_type"] == "tppl"
        self.df_trees.loc[tppl_idx, "tree_samples"] = [
            self.get_tppl_tree(
                row["filename"],
                burnin=burnin,
                subsample=subsample,
                end=end,
            )
            for _, row in self.df_fn.loc[tppl_idx, :].iterrows()
        ]
        self.active_idx = self.active_idx.union(self.tppl_idx)

    def get_tppl_tree(self, fn, burnin=0, subsample=1, end=None):
        tree = RepertoireTree()
        with open(fn) as fh:
            sliced_lines = islice(fh, burnin, end, subsample)
            for line in sliced_lines:
                sample = json.loads(line)
                if "__data__" in sample:
                    tree.append_json(sample["__data__"]["tree"])
        return tree

    def avg_trees(self, size=3, clean=True):
        def rep_to_vec(vec, size):
            nvec = np.array(vec)
            return np.eye(size)[nvec]

        def c_avg_rep(rep_arr, size=3):
            avg_rep = np.array([rep_to_vec(r, size) for r in rep_arr]).mean(axis=0)
            return avg_rep

        self.df_trees.loc[self.active_idx, "avg_tree"] = self.df_trees.loc[
            self.active_idx, "tree_samples"
        ].apply(
            lambda t: t.map(
                c_avg_rep,
                ["repertoire"],
                "avg_rep",
            )
        )
        if clean:
            self.df_trees = self.df_trees.drop(columns="tree_samples")

    def row_to_label(
        self,
        row,
        keys=[
            "genid",
            "param_id",
            "file_type",
            "model_dir",
            "model_name",
            "compile_id",
        ],
    ):
        return ".".join(str(row[key]) for key in keys)

    def heat_trees(self, imshow_kwargs={"cmap": "BuGn"}, rerender=False, clean=False):
        def mat_to_heatmap_file(mat, node_name, tempdir):
            tempdir.mkdir(parents=True, exist_ok=True)
            filename = f"avg_repertoire.{node_name}.png"
            file_path = tempdir / filename
            if not file_path.exists() or rerender:
                plt.figure(figsize=(2, 2), dpi=100)  # Fixed size and resolution
                plt.imshow(mat.T, aspect="auto", **imshow_kwargs)
                # plt.axis("off")  # Hide axes
                plt.xlabel("Host")
                plt.ylabel("State")
                plt.savefig(file_path, bbox_inches="tight", pad_inches=0)
                plt.close()

            return os.path.abspath(str(file_path))

        self.df_trees.loc[self.active_idx, "heat_map_tree"] = self.df_trees.loc[
            self.active_idx,
            :,
        ].apply(
            lambda row: row["avg_tree"].map(
                (
                    lambda mat, nn: mat_to_heatmap_file(
                        mat, nn, self.global_tempdir / self.row_to_label(row)
                    )
                ),
                ["avg_rep", "label"],
                "hmapfn",
            ),
            axis=1,
        )
        if clean:
            self.df_trees = self.df_trees.drop(columns="avg_tree")

    def dot_trees(self, clean=False):
        self.df_trees.loc[self.active_idx, "dot_tree"] = self.df_trees.loc[
            self.active_idx, "heat_map_tree"
        ].apply(lambda t: t.create_graphviz_tree())
        if clean:
            self.df_trees = self.df_trees.drop(columns="heat_map_tree")

    def render(self, rerender=False, debug=True):
        def render_tree(dot, tempdir, fn="tree.png", rerender=False):
            if debug:
                print(f"Rendering tree: {tempdir}")
            path = tempdir / fn
            if not path.exists() or rerender:
                path = dot.render(outfile=path)
            return path

        self.df_trees.loc[self.active_idx, "png_tree"] = self.df_trees.loc[
            self.active_idx,
        ].apply(
            lambda row: render_tree(
                row["dot_tree"],
                self.global_tempdir / self.row_to_label(row),
                rerender=rerender,
            ),
            axis=1,
        )

    def get_tree_images(self, render_tree=True, rerender=True):
        if "render_trees" not in self.df_trees.columns or rerender:
            self.avg_trees(clean=True)
            self.heat_trees(rerender=rerender, clean=True)
            self.dot_trees(clean=True)
            if render_tree:
                self.render(rerender=rerender)
        return self.df_trees
