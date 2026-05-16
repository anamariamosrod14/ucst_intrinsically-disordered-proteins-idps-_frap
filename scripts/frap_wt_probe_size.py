
# ============================================================
# FRAP ANALYSIS CORE
# Regions expected in Excel sheets:
#   Region 1 = pore ROI
#   Region 2 = matrix / continuous phase ROI
#   Region 3 = unbleached reference ROI
#
# Methods:
#   Double normalization:
#       Fdn_j(t) = (I_j(t)/<I_j,pre>) / (I_ref(t)/<I_ref,pre>)
#   Recovery normalization:
#       R_j(t) = (Fdn_j(t) - F0_j)/(Finf_j - F0_j)
#   Mobile fraction:
#       M_j = (Finf_j - F0_j)/(1 - F0_j)
#   Apparent diffusion:
#       D_app,j = gamma*w^2/(4*t1/2,j), gamma = 0.88
#   Relative diffusivity:
#       K_j = D_app,j/D_PBS
#
# Notes:
#   K is a relative diffusivity, not a hindrance factor that increases with restriction.
#   Lower K means stronger diffusive restriction relative to PBS.
# ============================================================

import os
import re
import math
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import kruskal, mannwhitneyu, rankdata, norm, wilcoxon

# -----------------------------
# USER SETTINGS
# -----------------------------
SEQ_ORDER = ["WT", "WS", "FS", "AS"]
SEQ_COLORS_BASE = {"WT": "#CC6677", "WS": "#4477AA", "FS": "#DDAA33", "AS": "#117733"}
SEQ_LABELS = {"WT": "GRGDSPYS", "WS": "GRGNSPWS", "FS": "GRGNSPFS", "AS": "GRGASPYA"}

# Filename tokens used to infer sequence identity
FNAME_TO_SEQ = {"GRGDSPYS": "WT", "GRGNSPWS": "WS", "GRGNSPFS": "FS", "GRGASPYA": "AS",
                "WT": "WT", "YS": "WT", "WS": "WS", "FS": "FS", "AS": "AS", "YA": "AS"}

# Probe sizes and official PBS diffusion values used consistently in this analysis
PROBES_TO_USE = [40, 250, 500]
D_FREE = {
    40:  {"D": 9.56, "D_sd": 0.62},
    250: {"D": 7.54, "D_sd": 0.61},
    500: {"D": 5.70, "D_sd": 1.98},
}

# FRAP parameters
ROI_DIAM_UM = 14.0
W_UM = ROI_DIAM_UM / 2.0
GAMMA = 0.88
N_LAST = 10
DT_GRID = 0.5
TMIN = 0.0

# Plot settings
FIG_W = 6.8
FIG_H = 4.8
LINE_W = 2.6
FONT_SIZE = 14
AXIS_W = 2.2
TICK_W = 1.8
TICK_LEN = 5
ALPHA_BAND = 0.18
BOX_WHIS = (0, 100)

# For publication, do not exclude high K values unless a documented QC reason exists.
APPLY_QC_FOR_K_ONLY = False
K_QC_MAX = 0.9

def lighten(color, amount=0.45):
    c = mcolors.to_rgb(color)
    return tuple(1 - (1 - x) * (1 - amount) for x in c)

SEQ_COLORS_PORE = {s: SEQ_COLORS_BASE[s] for s in SEQ_COLORS_BASE}
SEQ_COLORS_MATRIX = {s: lighten(SEQ_COLORS_BASE[s], amount=0.45) for s in SEQ_COLORS_BASE}

# -----------------------------
# General helpers
# -----------------------------
def _norm_name(s):
    return re.sub(r"\s+", " ", str(s).strip())

def find_col(df, candidates):
    cols = [_norm_name(c) for c in df.columns]
    for cand in candidates:
        if cand.startswith("re:"):
            pat = re.compile(cand[3:], re.IGNORECASE)
            for c in cols:
                if pat.search(c):
                    return c
        else:
            cand_n = _norm_name(cand)
            for c in cols:
                if c.lower() == cand_n.lower():
                    return c
    return None

def infer_probe_kDa(sheet_name):
    m = re.search(r"(\d+)\s*_?\s*kda", str(sheet_name), re.IGNORECASE)
    return int(m.group(1)) if m else np.nan

def infer_replicate(sheet_name):
    m = re.search(r"_(\d+)\s*$", str(sheet_name))
    return int(m.group(1)) if m else np.nan

def seq_from_filename(fname):
    key = os.path.basename(str(fname)).upper()
    for token, seq in FNAME_TO_SEQ.items():
        if token.upper() in key:
            return seq
    return "UNK"

# -----------------------------
# FRAP calculations
# -----------------------------
def double_normalization(t_rel, I_roi, I_ref):
    t_rel = np.asarray(t_rel, float)
    I_roi = np.asarray(I_roi, float)
    I_ref = np.asarray(I_ref, float)

    pre = t_rel < 0
    if np.sum(pre) < 3:
        nan = np.full_like(I_roi, np.nan, float)
        return nan, np.nan, np.nan

    roi_pre = np.nanmean(I_roi[pre])
    ref_pre = np.nanmean(I_ref[pre])
    if (not np.isfinite(roi_pre)) or (not np.isfinite(ref_pre)) or roi_pre == 0 or ref_pre == 0:
        nan = np.full_like(I_roi, np.nan, float)
        return nan, roi_pre, ref_pre

    Fdn = (I_roi / roi_pre) / (I_ref / ref_pre)
    return Fdn, roi_pre, ref_pre

def compute_R(t_rel, Fdn, n_last=N_LAST):
    t_rel = np.asarray(t_rel, float)
    Fdn = np.asarray(Fdn, float)
    if np.all(~np.isfinite(Fdn)):
        return np.full_like(Fdn, np.nan), np.nan, np.nan

    idx0 = int(np.argmin(np.abs(t_rel)))
    F0 = float(Fdn[idx0])
    if len(Fdn) < n_last:
        return np.full_like(Fdn, np.nan), F0, np.nan
    Finf = float(np.nanmean(Fdn[-n_last:]))

    if (not np.isfinite(F0)) or (not np.isfinite(Finf)) or np.isclose(Finf, F0):
        return np.full_like(Fdn, np.nan), F0, Finf

    R = (Fdn - F0) / (Finf - F0)
    return R, F0, Finf

def compute_mobile_fraction(F0, Finf):
    if (not np.isfinite(F0)) or (not np.isfinite(Finf)) or np.isclose(1.0, F0):
        return np.nan
    return (Finf - F0) / (1.0 - F0)

def compute_t_half(t_rel, R):
    t_rel = np.asarray(t_rel, float)
    R = np.asarray(R, float)
    m = np.isfinite(t_rel) & np.isfinite(R)
    t = t_rel[m]
    r = R[m]
    if len(t) < 2:
        return np.nan

    post = t >= 0
    t = t[post]
    r = r[post]
    if len(t) < 2:
        return np.nan

    idx = np.argsort(t)
    t = t[idx]
    r = r[idx]

    for i in range(1, len(r)):
        if r[i-1] < 0.5 <= r[i]:
            t1, t2 = t[i-1], t[i]
            r1, r2 = r[i-1], r[i]
            if np.isclose(r2, r1):
                return float(t2)
            return float(t1 + (0.5 - r1) * (t2 - t1) / (r2 - r1))
    return np.nan

def compute_D_app_from_t_half(t_half_s, w_um=W_UM, gamma=GAMMA):
    if (not np.isfinite(t_half_s)) or t_half_s <= 0:
        return np.nan
    return float(gamma * (w_um ** 2) / (4.0 * t_half_s))

def compute_K(probe_kDa, D_app):
    if not np.isfinite(probe_kDa) or not np.isfinite(D_app) or D_app <= 0:
        return np.nan
    probe_kDa = int(probe_kDa)
    if probe_kDa not in D_FREE:
        return np.nan
    return float(D_app) / float(D_FREE[probe_kDa]["D"])

# -----------------------------
# Statistics
# -----------------------------
def holm_adjust(pvals):
    pvals = np.asarray(pvals, float)
    m = len(pvals)
    out = np.full(m, np.nan, float)
    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out
    pv = pvals[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adj_ranked = np.empty(len(pv), float)
    running_max = 0.0
    for i, p in enumerate(ranked):
        val = (len(pv) - i) * p
        running_max = max(running_max, val)
        adj_ranked[i] = min(running_max, 1.0)
    adj_valid = np.empty(len(pv), float)
    adj_valid[order] = adj_ranked
    out[valid] = adj_valid
    return out

def dunn_posthoc(df, value_col, group_col):
    d = df[[group_col, value_col]].dropna().copy()
    d = d[np.isfinite(d[value_col])]
    groups = list(d[group_col].unique())
    if len(groups) < 2:
        return pd.DataFrame(columns=["group1", "group2", "z", "p_raw", "n1", "n2"])

    values = d[value_col].values
    ranks = rankdata(values, method="average")
    d["_rank"] = ranks
    N = len(d)

    _, counts = np.unique(values, return_counts=True)
    T = np.sum(counts ** 3 - counts)
    C = 1.0 - T / (N ** 3 - N) if N > 1 else 1.0
    if C <= 0:
        C = 1.0

    grp = d.groupby(group_col)["_rank"].agg(["count", "sum"]).rename(columns={"count": "n", "sum": "R"})
    rows = []
    for g1, g2 in itertools.combinations(groups, 2):
        n1, R1 = grp.loc[g1, "n"], grp.loc[g1, "R"]
        n2, R2 = grp.loc[g2, "n"], grp.loc[g2, "R"]
        mean_r1, mean_r2 = R1 / n1, R2 / n2
        denom = np.sqrt((N * (N + 1) / 12.0) * C * (1.0 / n1 + 1.0 / n2))
        if denom == 0:
            z = np.nan
            p_raw = 1.0
        else:
            z = (mean_r1 - mean_r2) / denom
            p_raw = 2.0 * norm.sf(abs(z))
        rows.append(dict(group1=g1, group2=g2, z=z, p_raw=p_raw, n1=int(n1), n2=int(n2)))
    return pd.DataFrame(rows)

def kw_plus_dunn_holm(df, value_col, group_col, min_n_per_group=2):
    d = df[[group_col, value_col]].dropna().copy()
    d = d[np.isfinite(d[value_col])]
    sizes = d.groupby(group_col)[value_col].count()
    keep = sizes[sizes >= min_n_per_group].index.tolist()
    d = d[d[group_col].isin(keep)].copy()
    groups = sorted(d[group_col].unique())

    if len(groups) < 2:
        kw = pd.DataFrame({"H": [np.nan], "p": [np.nan], "k": [len(groups)]})
        pairs = pd.DataFrame(columns=["group1", "group2", "z", "p_raw", "p_adj", "n1", "n2"])
        return kw, pairs

    arrays = [d.loc[d[group_col] == g, value_col].values for g in groups]
    H, p_kw = kruskal(*arrays)
    kw = pd.DataFrame({"H": [float(H)], "p": [float(p_kw)], "k": [int(len(groups))]})

    if len(groups) == 2:
        g1, g2 = groups
        x = d.loc[d[group_col] == g1, value_col].values
        y = d.loc[d[group_col] == g2, value_col].values
        stat, p_raw = mannwhitneyu(x, y, alternative="two-sided")
        pairs = pd.DataFrame([{"group1": g1, "group2": g2, "z": np.nan, "p_raw": float(p_raw),
                               "p_adj": float(p_raw), "n1": int(len(x)), "n2": int(len(y))}])
        return kw, pairs

    pairs = dunn_posthoc(d, value_col=value_col, group_col=group_col)
    if not pairs.empty:
        pairs["p_adj"] = holm_adjust(pairs["p_raw"].values)
    else:
        pairs["p_adj"] = []
    return kw, pairs

def p_to_symbol(p):
    if not np.isfinite(p):
        return "na"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"

# -----------------------------
# Plot helpers
# -----------------------------
def stylize_axes(ax):
    ax.grid(False)
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(AXIS_W)
    ax.tick_params(width=TICK_W, length=TICK_LEN, labelsize=FONT_SIZE)

def build_common_grid(rows, ykey, dt=DT_GRID, tmin=TMIN, tmax=None):
    tmaxs = []
    for r in rows:
        t = np.asarray(r["t_rel"], float)
        y = np.asarray(r[ykey], float)
        m = np.isfinite(t) & np.isfinite(y) & (t >= 0)
        if np.sum(m) >= 2:
            tmaxs.append(np.max(t[m]))
    if not tmaxs:
        return None
    if tmax is None:
        tmax = float(np.min(tmaxs))
    if tmax <= tmin + dt:
        return None
    return np.arange(tmin, tmax + 1e-9, dt)

def interp_to_grid(t, y, grid):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(y)
    t = t[m]
    y = y[m]
    if len(t) < 2:
        return np.full_like(grid, np.nan, float)
    idx = np.argsort(t)
    return np.interp(grid, t[idx], y[idx], left=np.nan, right=np.nan)

def mean_sd_on_grid(rows, ykey, grid):
    Y = np.vstack([interp_to_grid(r["t_rel"], r[ykey], grid) for r in rows])
    return np.nanmean(Y, axis=0), np.nanstd(Y, axis=0)

def ylab(metric_col):
    if metric_col in ["t_half_s", "t_half"]:
        return r"$t_{1/2}$ (s)"
    if metric_col in ["D_app_um2_s", "D_app"]:
        return r"$D_{\mathrm{app}}$ (µm$^2$/s)"
    if metric_col in ["mobile_fraction", "M"]:
        return "Mobile fraction"
    if metric_col == "K":
        return r"$K = D_{\mathrm{app}}/D_{\mathrm{PBS}}$"
    return metric_col

def save_or_show(fig, path=None, show=True):
    if path:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def draw_pairwise_brackets(ax, pairs_df, group_order, data_arrays, symbol_col="symbol", show_ns=True):
    """Draw pairwise significance brackets on a boxplot.

    Parameters
    ----------
    ax : matplotlib axis
    pairs_df : DataFrame with group1, group2, and symbol/p_adj columns
    group_order : list of group identifiers in the x-axis order
    data_arrays : list of arrays corresponding to group_order
    symbol_col : column containing '*', '**', '***', or 'ns'
    show_ns : if True, draw ns labels as well as significant comparisons
    """
    if pairs_df is None or pairs_df.empty:
        return

    x_pos = {g: i + 1 for i, g in enumerate(group_order)}

    finite_vals = []
    for arr in data_arrays:
        arr = np.asarray(arr, dtype=float)
        finite_vals.extend(arr[np.isfinite(arr)].tolist())
    if len(finite_vals) == 0:
        return

    y_min = float(np.nanmin(finite_vals))
    y_max = float(np.nanmax(finite_vals))
    y_range = y_max - y_min
    if (not np.isfinite(y_range)) or y_range <= 0:
        y_range = max(abs(y_max), 1.0)

    # Draw shorter comparisons first so nested brackets stack cleanly.
    tmp = pairs_df.copy()
    tmp = tmp[tmp["group1"].isin(group_order) & tmp["group2"].isin(group_order)].copy()
    if tmp.empty:
        return
    tmp["_span"] = tmp.apply(lambda r: abs(x_pos[r["group2"]] - x_pos[r["group1"]]), axis=1)
    tmp = tmp.sort_values(["_span", "p_adj" if "p_adj" in tmp.columns else "p_raw"], ascending=[True, True])

    y = y_max + 0.12 * y_range
    step = 0.16 * y_range
    cap = 0.04 * y_range

    for _, row in tmp.iterrows():
        sym = row.get(symbol_col, None)
        if sym is None or (isinstance(sym, float) and np.isnan(sym)):
            sym = p_to_symbol(row.get("p_adj", np.nan))
        if (sym == "ns") and (not show_ns):
            continue

        x1, x2 = x_pos[row["group1"]], x_pos[row["group2"]]
        if x1 > x2:
            x1, x2 = x2, x1
        ax.plot([x1, x1, x2, x2], [y, y + cap, y + cap, y], lw=1.6, c="black", clip_on=False)
        ax.text((x1 + x2) / 2, y + cap + 0.02 * y_range, str(sym),
                ha="center", va="bottom", fontsize=FONT_SIZE - 1, clip_on=False)
        y += step

    bottom, top = ax.get_ylim()
    if y + step > top:
        ax.set_ylim(bottom, y + step)

# ============================================================
# WT-ONLY PROBE SIZE PIPELINE
# Compares 40, 250, and 500 kDa probes in WT GRGDSPYS.
# ============================================================

PROBE_COLORS = {40: "#CC6677", 250: "#8A4653", 500: "#5A2D36"}

def analyze_wt_excel_file(fname):
    metrics_rows = []
    point_tables = []
    runs = []

    xls = pd.ExcelFile(fname)
    sheet_counter = 0
    for sheet in xls.sheet_names:
        probe = infer_probe_kDa(sheet)
        if not np.isfinite(probe) or int(probe) not in PROBES_TO_USE:
            continue
        sheet_counter += 1
        rep = infer_replicate(sheet)
        if not np.isfinite(rep):
            rep = sheet_counter

        df = pd.read_excel(fname, sheet_name=sheet)
        df.columns = [_norm_name(c) for c in df.columns]

        col_trel = find_col(df, ["Adjusted time [s]", "re:adjusted.*time"])
        col_r1 = find_col(df, ["Intensity Region 1", "re:intensity\\s*region\\s*1", "re:region\\s*1"])
        col_r2 = find_col(df, ["Intensity Region 2", "re:intensity\\s*region\\s*2", "re:region\\s*2"])
        col_r3 = find_col(df, ["Intensity Region 3", "re:intensity\\s*region\\s*3", "re:region\\s*3"])

        if any(c is None for c in [col_trel, col_r1, col_r2, col_r3]):
            print(f"[SKIP] Missing required columns in sheet {sheet}")
            continue

        t_rel = df[col_trel].to_numpy(float)
        I_pore = df[col_r1].to_numpy(float)
        I_matrix = df[col_r2].to_numpy(float)
        I_ref = df[col_r3].to_numpy(float)

        Fdn_pore, _, _ = double_normalization(t_rel, I_pore, I_ref)
        Fdn_matrix, _, _ = double_normalization(t_rel, I_matrix, I_ref)

        R_pore, F0_pore, Finf_pore = compute_R(t_rel, Fdn_pore)
        R_matrix, F0_matrix, Finf_matrix = compute_R(t_rel, Fdn_matrix)

        t_half_pore = compute_t_half(t_rel, R_pore)
        t_half_matrix = compute_t_half(t_rel, R_matrix)
        D_pore = compute_D_app_from_t_half(t_half_pore)
        D_matrix = compute_D_app_from_t_half(t_half_matrix)
        M_pore = compute_mobile_fraction(F0_pore, Finf_pore)
        M_matrix = compute_mobile_fraction(F0_matrix, Finf_matrix)

        runs.append(dict(sheet=sheet, probe_kDa=probe, replicate=rep, t_rel=t_rel,
                         R_pore=R_pore, R_matrix=R_matrix,
                         Fdn_pore=Fdn_pore, Fdn_matrix=Fdn_matrix))

        for roi_type, t_half, D_app, M in [
            ("pore", t_half_pore, D_pore, M_pore),
            ("matrix", t_half_matrix, D_matrix, M_matrix),
        ]:
            metrics_rows.append(dict(sequence="WT", sheet=sheet, probe_kDa=int(probe), replicate=rep,
                                     roi_type=roi_type, t_half_s=t_half, D_app_um2_s=D_app,
                                     mobile_fraction=M, K=compute_K(probe, D_app)))

        point_tables.append(pd.DataFrame({
            "sheet": sheet,
            "probe_kDa": int(probe),
            "replicate": rep,
            "Time_rel_s": t_rel,
            "Fdn_pore": Fdn_pore,
            "Fdn_matrix": Fdn_matrix,
            "R_pore": R_pore,
            "R_matrix": R_matrix
        }))

    metrics_df = pd.DataFrame(metrics_rows)
    points_df = pd.concat(point_tables, ignore_index=True) if point_tables else pd.DataFrame()
    return metrics_df, points_df, runs

def plot_wt_recovery(runs, roi_type, output_dir=None, show=True):
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ykey = "R_pore" if roi_type == "pore" else "R_matrix"
    ls = "-" if roi_type == "pore" else "--"

    for probe in PROBES_TO_USE:
        rows = [{"t_rel": r["t_rel"], ykey: r[ykey]} for r in runs if int(r["probe_kDa"]) == int(probe)]
        if not rows:
            continue
        grid = build_common_grid(rows, ykey)
        if grid is None:
            continue
        mu, sd = mean_sd_on_grid(rows, ykey, grid)
        ax.plot(grid, mu, color=PROBE_COLORS.get(probe, "#333333"), lw=LINE_W, ls=ls, label=f"{probe} kDa")
        ax.fill_between(grid, mu - sd, mu + sd, color=PROBE_COLORS.get(probe, "#333333"), alpha=ALPHA_BAND)

    ax.set_xlim(left=0)
    ax.set_xlabel("Time (s)", fontsize=FONT_SIZE)
    ax.set_ylabel("Normalized Fluorescence Recovery", fontsize=FONT_SIZE)
    ax.set_title(f"WT GRGDSPYS — {roi_type} recovery", fontsize=FONT_SIZE)
    ax.legend(frameon=False, loc="lower right", fontsize=FONT_SIZE-2)
    stylize_axes(ax)
    plt.tight_layout()

    path = None
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f"WT_recovery_{roi_type}.png")
    save_or_show(fig, path=path, show=show)

def boxplot_probe_metric(metrics_df, metric_col, roi_type, output_dir=None, show=True):
    sub = metrics_df[(metrics_df["roi_type"] == roi_type) & np.isfinite(metrics_df[metric_col])].copy()
    probes_present = [p for p in PROBES_TO_USE if p in sub["probe_kDa"].unique()]
    if len(probes_present) < 2:
        return None, None

    data = [sub[sub["probe_kDa"] == p][metric_col].dropna().values for p in probes_present]
    labels = [f"{p} kDa" for p in probes_present]
    colors = [PROBE_COLORS.get(p, "#777777") for p in probes_present]

    kw, pairs = kw_plus_dunn_holm(sub.rename(columns={metric_col: "value"}), "value", "probe_kDa")
    if not pairs.empty:
        pairs["symbol"] = pairs["p_adj"].apply(p_to_symbol)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, whis=BOX_WHIS)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
        patch.set_linewidth(1.8)
    for key in ["whiskers", "caps", "medians"]:
        for line in bp[key]:
            line.set_linewidth(1.6)

    rng = np.random.default_rng(123)
    for i, vals in enumerate(data, 1):
        ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=38,
                   facecolors="white", edgecolors="black", linewidths=1.0, zorder=3)

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=FONT_SIZE)
    ax.set_ylabel(ylab(metric_col), fontsize=FONT_SIZE)
    ax.set_title(f"WT GRGDSPYS — {roi_type} — {ylab(metric_col)}", fontsize=FONT_SIZE)
    ax.set_ylim(bottom=0)
    draw_pairwise_brackets(ax, pairs, probes_present, data, symbol_col="symbol", show_ns=True)
    stylize_axes(ax)
    plt.tight_layout()

    path = None
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f"WT_box_{metric_col}_{roi_type}.png")
    save_or_show(fig, path=path, show=show)

    return kw, pairs

def run_wt_probe_size_pipeline(fname, output_dir="WT_FRAP_probe_size_output", show=True):
    os.makedirs(output_dir, exist_ok=True)
    metrics_df, points_df, runs = analyze_wt_excel_file(fname)

    if metrics_df.empty:
        raise RuntimeError("No WT FRAP metrics were extracted. Check sheet names and column names.")

    summary = metrics_df.groupby(["probe_kDa", "roi_type"]).agg(
        t_half_mean=("t_half_s", "mean"), t_half_sd=("t_half_s", "std"),
        D_app_mean=("D_app_um2_s", "mean"), D_app_sd=("D_app_um2_s", "std"),
        K_mean=("K", "mean"), K_sd=("K", "std"),
        M_mean=("mobile_fraction", "mean"), M_sd=("mobile_fraction", "std"),
        n=("t_half_s", "count")
    ).reset_index()

    kw_rows = []
    pairwise_rows = []
    for roi in ["pore", "matrix"]:
        plot_wt_recovery(runs, roi, output_dir=os.path.join(output_dir, "figures"), show=show)

    for metric in ["t_half_s", "D_app_um2_s", "K", "mobile_fraction"]:
        for roi in ["pore", "matrix"]:
            kw, pairs = boxplot_probe_metric(metrics_df, metric, roi, output_dir=os.path.join(output_dir, "figures"), show=show)
            if kw is not None:
                kw_rows.append(dict(metric=metric, roi_type=roi, H=kw["H"].iloc[0],
                                    p_global=kw["p"].iloc[0], k_groups=kw["k"].iloc[0]))
            if pairs is not None and not pairs.empty:
                for _, r in pairs.iterrows():
                    pairwise_rows.append(dict(metric=metric, roi_type=roi,
                                              group1=r["group1"], group2=r["group2"],
                                              z=r.get("z", np.nan), p_raw=r["p_raw"], p_adj=r["p_adj"],
                                              symbol=r.get("symbol", p_to_symbol(r["p_adj"])),
                                              n1=r.get("n1", np.nan), n2=r.get("n2", np.nan)))

    kw_df = pd.DataFrame(kw_rows)
    pairwise_df = pd.DataFrame(pairwise_rows)

    out_xlsx = os.path.join(output_dir, "WT_FRAP_probe_size_stats.xlsx")
    with pd.ExcelWriter(out_xlsx) as writer:
        metrics_df.to_excel(writer, index=False, sheet_name="metrics_per_rep")
        points_df.to_excel(writer, index=False, sheet_name="all_timepoints_long")
        summary.to_excel(writer, index=False, sheet_name="summary_means_sd")
        kw_df.to_excel(writer, index=False, sheet_name="KW_global")
        pairwise_df.to_excel(writer, index=False, sheet_name="pairwise_DunnHolm")

    print(f"Saved output folder: {output_dir}")
    print(f"Saved Excel summary: {out_xlsx}")
    return metrics_df, points_df, summary, kw_df, pairwise_df

# ============================================================
# GITHUB / LOCAL ENTRY POINT
# Usage:
#   python scripts/frap_wt_probe_size.py --excel_file data/raw/WT_FRAP_raw.xlsx --output_dir results/wt_probe_size
# ============================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="WT GRGDSPYS FRAP probe-size analysis for 40, 250, and 500 kDa probes."
    )
    parser.add_argument(
        "--excel_file",
        required=True,
        help="Excel file containing WT GRGDSPYS FRAP sheets. Sheet names should include probe size, e.g. YS_40kDa_1."
    )
    parser.add_argument(
        "--output_dir",
        default="WT_FRAP_probe_size_output",
        help="Directory where figures, tables, and Excel outputs will be saved."
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots interactively. By default, plots are saved and closed."
    )

    args = parser.parse_args()
    run_wt_probe_size_pipeline(
        fname=args.excel_file,
        output_dir=args.output_dir,
        show=args.show
    )
