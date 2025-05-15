import re
import pandas as pd
from datetime import timedelta


def parse_trace(tracefile):
    df_trace = pd.read_csv(tracefile, sep="\t")
    # Identify rows where status is not null – these are your main data rows
    df_trace["is_main"] = df_trace["status"].notnull()

    # Forward-fill a group identifier from the main rows
    df_trace["group"] = df_trace["task_id"].where(df_trace["is_main"]).ffill()

    # Separate main rows and note rows
    main_rows = df_trace[df_trace["is_main"]].copy()
    note_rows = df_trace[~df_trace["is_main"]].copy()

    # Aggregate notes by group
    notes_by_group = (
        note_rows.groupby("group")["task_id"]
        .apply(lambda x: "".join(list(map(lambda s: s.strip(), x))))
        .to_dict()
    )

    # Attach notes to the main rows
    main_rows["notes"] = main_rows["task_id"].map(notes_by_group)
    run_rows = main_rows[main_rows["name"].astype(str).str.startswith("run_hostrep")]
    run_rows["parsed_time"] = run_rows["realtime"].apply(parse_time)

    run_rows[["genid", "param_id", "compile_id", "file_type"]] = run_rows[
        "notes"
    ].apply(parse_model)
    df_trace = run_rows[
        [
            "genid",
            "param_id",
            "file_type",
            "compile_id",
            "status",
            "parsed_time",
            "%cpu",
        ]
    ]
    return df_trace


def parse_time(time_str):
    # Extract hours, minutes, seconds if present
    pattern = r"(?:(\d+)h)?\s*(?:(\d+)m)?\s*(?:(\d+)s)?"
    match = re.match(pattern, time_str.strip())
    if not match or time_str == "-":
        return None

    hours = int(match.group(1)) if match.group(1) else 0
    minutes = int(match.group(2)) if match.group(2) else 0
    seconds = int(match.group(3)) if match.group(3) else 0

    return timedelta(hours=hours, minutes=minutes, seconds=seconds).total_seconds()


def parse_model(script_str):
    match = re.search(r"output\.(\d+)\.(\d+)\.(\d+)\.json", script_str)
    if not match:
        file_type = "rb"
        match = re.search(r"out\.(\d+)\.(\d+)\.(\d+)", script_str)
    else:
        file_type = "tppl"
    if not match:
        return pd.Series([None, None, None, None])
    param_id = int(match.group(1))
    genid = int(match.group(2))
    compile_id = int(match.group(3))
    return pd.Series([genid, param_id, compile_id, file_type])
