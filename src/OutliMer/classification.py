from __future__ import annotations

import argparse
import html
import inspect
import json
import math
import os
import sys
from datetime import datetime, timezone
from typing import Any, Iterable, Sequence


def _load_pandas_numpy():
    try:
        import numpy as np
        import pandas as pd
    except Exception as exc:
        raise RuntimeError(
            "OutliMer classification requires pandas and numpy. "
            "Install OutliMer with its runtime dependencies."
        ) from exc
    return pd, np


def _load_plotting():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception as exc:
        raise RuntimeError(
            "Plotting requires matplotlib and seaborn. "
            "Use --no-plots to skip plot generation."
        ) from exc
    return plt, sns


def _load_sklearn():
    try:
        from sklearn.decomposition import PCA
        from sklearn.ensemble import IsolationForest
        from sklearn.neighbors import LocalOutlierFactor
    except Exception as exc:
        raise RuntimeError(
            "Union-mode anomaly scoring requires scikit-learn."
        ) from exc
    return PCA, IsolationForest, LocalOutlierFactor


def _load_sklearn_manifold():
    try:
        from sklearn.manifold import MDS
    except Exception as exc:
        raise RuntimeError("Jaccard MDS plotting requires scikit-learn.") from exc
    return MDS


def _load_scipy_clustering():
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist
    except Exception as exc:
        raise RuntimeError(
            "Dendrogram generation requires scipy. Use --no-plots to skip it."
        ) from exc
    return pdist, linkage, dendrogram


def _load_fisher_exact():
    try:
        from scipy.stats import fisher_exact
    except Exception as exc:
        raise RuntimeError("Hash enrichment requires scipy.") from exc
    return fisher_exact


def _validate_contamination(value: float) -> float:
    if value <= 0 or value > 0.5:
        raise ValueError("--contamination must be > 0 and <= 0.5")
    return value


def _parse_sample_list(value: str | None) -> list[str] | None:
    if not value:
        return None
    samples = [part.strip() for part in value.split(",") if part.strip()]
    return samples or None


def load_sample_metadata(path: str, sample_column: str = "sample"):
    """Load sample metadata from CSV/TSV keyed by sample name."""
    pd, _ = _load_pandas_numpy()
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    delimiter = "\t" if path.lower().endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=delimiter, dtype=str).fillna("")
    if sample_column not in df.columns:
        raise ValueError(f"metadata sample column {sample_column!r} not found")
    df = df.drop_duplicates(subset=[sample_column], keep="last")
    df = df.set_index(sample_column)
    return df


def attach_metadata(df: Any, metadata: Any | None) -> Any:
    if metadata is None or df.empty:
        return df
    joined = df.join(metadata, how="left")
    return joined


def select_samples_by_query(metadata: Any, query: str) -> list[str]:
    """Select samples from metadata using pandas query syntax."""
    if metadata is None:
        raise ValueError("metadata is required for metadata query selection")
    try:
        selected = metadata.query(query, engine="python")
    except Exception as exc:
        raise ValueError(f"metadata query failed: {query!r}: {exc}") from exc
    return [str(idx) for idx in selected.index]


def load_union_csv(path: str) -> Any:
    """Load a top-union summary CSV as rows=hashes, columns=samples."""
    pd, _ = _load_pandas_numpy()
    if not path or not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, index_col=0)
    for column in ("sequence", "total_count"):
        if column in df.columns:
            df = df.drop(columns=[column])
    return df.fillna(0).apply(pd.to_numeric, errors="coerce").fillna(0)


def load_report_csv(path: str) -> Any:
    """Load the per-sample report written by the main OutliMer CLI."""
    pd, _ = _load_pandas_numpy()
    if not path or not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, index_col=0)
    df.columns = [str(c).lower() for c in df.columns]
    expected = {"n_hashes", "n_in_db", "pct_in_db", "n_new_hashes"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"Report CSV missing expected columns: {sorted(missing)} "
            f"(found {sorted(df.columns)})"
        )
    for column in expected:
        df[column] = df[column].astype(float)
    return df


def explain_report_based(report_df: Any) -> Any:
    """Rank samples from the report CSV using simple overlap heuristics."""
    pd, _ = _load_pandas_numpy()
    df = report_df.copy()
    pct_comp = 1.0 - df["pct_in_db"]
    max_new = max(1.0, float(df["n_new_hashes"].max()))
    new_comp = df["n_new_hashes"] / max_new
    max_size = max(1.0, float(df["n_hashes"].max()))
    size_comp = 1.0 - (df["n_hashes"] / max_size)

    raw_score = (0.6 * pct_comp) + (0.3 * new_comp) + (0.1 * size_comp)
    span = float(raw_score.max() - raw_score.min())
    if span <= 1e-12:
        score = raw_score * 0.0
    else:
        score = (raw_score - raw_score.min()) / span

    reasons: list[str] = []
    median_new = float(df["n_new_hashes"].median())
    for _, row in df.iterrows():
        parts = []
        if row["pct_in_db"] < 0.2:
            parts.append(f"low overlap with DB (pct_in_db={row['pct_in_db']:.3f})")
        elif row["pct_in_db"] < 0.5:
            parts.append(
                f"moderate overlap with DB (pct_in_db={row['pct_in_db']:.3f})"
            )
        else:
            parts.append(f"high overlap with DB (pct_in_db={row['pct_in_db']:.3f})")

        if row["n_new_hashes"] > median_new * 2:
            parts.append(f"many novel hashes (n_new_hashes={int(row['n_new_hashes'])})")
        elif row["n_new_hashes"] > 0:
            parts.append(f"some novel hashes (n_new_hashes={int(row['n_new_hashes'])})")

        if row["n_hashes"] < max(50, int(max_size * 0.01)):
            parts.append(
                f"very small sample (n_hashes={int(row['n_hashes'])}); may be noisy"
            )
        reasons.append("; ".join(parts))

    out = pd.DataFrame({
        "sample": df.index,
        "pct_in_db": df["pct_in_db"].values,
        "n_hashes": df["n_hashes"].values,
        "n_new_hashes": df["n_new_hashes"].values,
        "reason_score": score.values,
        "reasons": reasons,
    }).set_index("sample")
    return out.sort_values("reason_score", ascending=False)


def prepare_feature_matrix(df: Any, top_M: int = 2000, min_samples: int = 1):
    """Return raw, binary, and log1p sample-by-hash feature matrices."""
    _, np = _load_pandas_numpy()
    if min_samples < 1:
        raise ValueError("--min-samples must be >= 1")
    if top_M < 0:
        raise ValueError("--top-features must be >= 0")

    X = df.T.astype(float).fillna(0)
    prevalence = X.gt(0).sum(axis=0)
    X = X.loc[:, prevalence >= min_samples]
    if top_M and X.shape[1] > top_M:
        totals = X.sum(axis=0)
        X = X.loc[:, totals.sort_values(ascending=False).head(top_M).index]
    X_binary = (X > 0).astype(int)
    X_log = np.log1p(X)
    return X, X_binary, X_log


def _ranked_label_samples(
    ranked: Any | None,
    max_labels: int = 12,
    extra_samples: Iterable[str] | None = None,
) -> set[str]:
    if max_labels <= 0:
        return set()
    labels = set(str(s) for s in (extra_samples or []))
    if ranked is not None and getattr(ranked, "shape", (0,))[0]:
        if "anomaly_score" in ranked.columns:
            ranked_samples = ranked.sort_values(
                "anomaly_score",
                ascending=False,
            ).head(max_labels).index
        else:
            ranked_samples = ranked.head(max_labels).index
        labels.update(str(s) for s in ranked_samples)
    return labels


def _offset_for_label(x_values: Any, y_values: Any) -> tuple[float, float]:
    _, np = _load_pandas_numpy()
    x_span = float(np.nanmax(x_values) - np.nanmin(x_values)) if len(x_values) else 1.0
    y_span = float(np.nanmax(y_values) - np.nanmin(y_values)) if len(y_values) else 1.0
    return (x_span or 1.0) * 0.01, (y_span or 1.0) * 0.01


def compute_pca_plot(
    X_log: Any,
    out_path: str,
    ranked: Any | None = None,
    max_labels: int = 12,
) -> None:
    n_samples, n_features = X_log.shape
    if n_samples < 2 or n_features < 2:
        print(
            "Skipping PCA plot: need >=2 samples and >=2 features "
            f"(have {n_samples} samples, {n_features} features)"
        )
        return
    PCA, _, _ = _load_sklearn()
    plt, sns = _load_plotting()
    n_comp = min(4, n_samples, n_features)
    pca = PCA(n_components=n_comp)
    Z = pca.fit_transform(X_log.values)
    fig, ax = plt.subplots(figsize=(7, 6))
    scores = None
    if ranked is not None and "anomaly_score" in ranked.columns:
        scores = ranked.reindex(X_log.index)["anomaly_score"].astype(float).fillna(0)
    if scores is not None:
        points = ax.scatter(
            Z[:, 0],
            Z[:, 1],
            c=scores.values,
            s=70,
            cmap="viridis",
            edgecolors="black",
            linewidths=0.3,
            alpha=0.85,
        )
        fig.colorbar(points, ax=ax, label="Anomaly score")
    else:
        sns.scatterplot(x=Z[:, 0], y=Z[:, 1], s=60, ax=ax)
    labels = _ranked_label_samples(ranked, max_labels=max_labels)
    if not labels and n_samples <= max_labels:
        labels = set(str(s) for s in X_log.index)
    dx, dy = _offset_for_label(Z[:, 0], Z[:, 1])
    for i, sample in enumerate(X_log.index):
        if str(sample) in labels:
            ax.text(Z[i, 0] + dx, Z[i, 1] + dy, sample, fontsize=8)
    variance = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1] * 100:.1f}%)")
    ax.set_title("PCA (log1p counts)")
    ax.margins(x=0.08, y=0.08)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_jaccard_distance_matrix(X_binary: Any) -> Any:
    """Return a square sample-by-sample Jaccard distance matrix."""
    pd, np = _load_pandas_numpy()
    arr = X_binary.values.astype(bool)
    n_samples = arr.shape[0]
    distances = np.zeros((n_samples, n_samples), dtype=float)
    if n_samples <= 1 or arr.shape[1] == 0:
        return pd.DataFrame(distances, index=X_binary.index, columns=X_binary.index)
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            union = int(np.logical_or(arr[i], arr[j]).sum())
            if union == 0:
                distance = 0.0
            else:
                intersection = int(np.logical_and(arr[i], arr[j]).sum())
                distance = 1.0 - (intersection / union)
            distances[i, j] = distances[j, i] = distance
    return pd.DataFrame(distances, index=X_binary.index, columns=X_binary.index)


def compute_mean_jaccard_distance(X_binary: Any) -> Any:
    """Compute mean pairwise Jaccard distance without requiring scipy."""
    pd, _ = _load_pandas_numpy()
    distances = compute_jaccard_distance_matrix(X_binary)
    if distances.empty:
        return pd.Series([], index=X_binary.index, dtype=float)
    return pd.Series(distances.values.mean(axis=1), index=X_binary.index)


def compute_sample_qc_metrics(df: Any, ranked: Any | None = None) -> Any:
    """Compute per-sample hash-depth and cohort-uniqueness diagnostics."""
    pd, np = _load_pandas_numpy()
    counts = df.fillna(0).astype(float)
    presence = counts.gt(0)
    prevalence = presence.sum(axis=1)
    n_hashes = presence.sum(axis=0).astype(int)
    unique_hashes = presence.loc[prevalence == 1].sum(axis=0).astype(int)
    shared_hashes = n_hashes - unique_hashes
    total_count = counts.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_unique = (unique_hashes / n_hashes.replace(0, np.nan)).fillna(0.0)
    qc = pd.DataFrame({
        "n_hashes": n_hashes,
        "total_count": total_count,
        "unique_hashes": unique_hashes,
        "shared_hashes": shared_hashes,
        "pct_unique": pct_unique,
        "pct_shared": 1.0 - pct_unique,
    })
    if ranked is not None and "anomaly_score" in ranked.columns:
        qc = qc.join(ranked[["anomaly_score"]], how="left")
        qc["anomaly_score"] = qc["anomaly_score"].fillna(0.0)
    return qc


def compute_isolation_forest(X_log: Any, contamination: float = 0.05) -> Any:
    pd, _ = _load_pandas_numpy()
    if X_log.shape[0] <= 1 or X_log.shape[1] == 0:
        return pd.Series([0.0] * X_log.shape[0], index=X_log.index, dtype=float)
    _validate_contamination(contamination)
    _, IsolationForest, _ = _load_sklearn()
    clf = IsolationForest(
        n_estimators=200,
        contamination=contamination,
        random_state=0,
    )
    clf.fit(X_log.values)
    return pd.Series(-clf.score_samples(X_log.values), index=X_log.index)


def compute_lof(X_log: Any, contamination: float = 0.05) -> Any:
    pd, _ = _load_pandas_numpy()
    if X_log.shape[0] <= 1 or X_log.shape[1] == 0:
        return pd.Series([0.0] * X_log.shape[0], index=X_log.index, dtype=float)
    _validate_contamination(contamination)
    _, _, LocalOutlierFactor = _load_sklearn()
    n_neighbors = max(1, min(20, X_log.shape[0] - 1))
    lof = LocalOutlierFactor(
        n_neighbors=n_neighbors,
        contamination=contamination,
        metric="euclidean",
    )
    lof.fit_predict(X_log.values)
    return pd.Series(-lof.negative_outlier_factor_, index=X_log.index)


def compute_dendrogram(X_binary: Any, out_path: str) -> None:
    if X_binary.shape[0] <= 1 or X_binary.shape[1] == 0:
        return
    pdist, linkage, dendrogram = _load_scipy_clustering()
    plt, _ = _load_plotting()
    dist_vec = pdist(X_binary.values, metric="jaccard")
    Z = linkage(dist_vec, method="average")
    fig = plt.figure(figsize=(10, 6))
    dendrogram(Z, labels=X_binary.index.tolist(), leaf_rotation=90)
    plt.title("Hierarchical clustering (Jaccard)")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_clustered_jaccard_heatmap(
    X_binary: Any,
    out_path: str,
) -> None:
    """Write a clustered sample-by-sample Jaccard similarity heatmap."""
    if X_binary.shape[0] <= 1 or X_binary.shape[1] == 0:
        return
    pdist, linkage, _ = _load_scipy_clustering()
    plt, sns = _load_plotting()
    dist_vec = pdist(X_binary.values, metric="jaccard")
    Z = linkage(dist_vec, method="average")
    distances = compute_jaccard_distance_matrix(X_binary)
    similarity = 1.0 - distances
    size = max(7.0, min(16.0, X_binary.shape[0] * 0.13))
    grid = sns.clustermap(
        similarity,
        row_linkage=Z,
        col_linkage=Z,
        cmap="viridis",
        vmin=0,
        vmax=1,
        figsize=(size, size),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Jaccard similarity"},
    )
    grid.fig.suptitle("Clustered Jaccard Similarity", y=1.02)
    grid.fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(grid.fig)


def compute_jaccard_mds_plot(
    X_binary: Any,
    out_path: str,
    ranked: Any | None = None,
    max_labels: int = 12,
) -> None:
    """Write a two-dimensional MDS plot from Jaccard distances."""
    if X_binary.shape[0] <= 2 or X_binary.shape[1] == 0:
        return
    _, np = _load_pandas_numpy()
    plt, _ = _load_plotting()
    MDS = _load_sklearn_manifold()
    distances = compute_jaccard_distance_matrix(X_binary)
    mds_params = inspect.signature(MDS).parameters
    mds_kwargs: dict[str, object] = {
        "n_components": 2,
        "random_state": 0,
    }
    if "normalized_stress" in mds_params:
        mds_kwargs["normalized_stress"] = "auto"
    if "init" in mds_params:
        mds_kwargs["init"] = "random"
    if "metric_mds" in mds_params:
        mds_kwargs["metric_mds"] = True
        mds_kwargs["metric"] = "precomputed"
    else:
        mds_kwargs["dissimilarity"] = "precomputed"
    mds = MDS(**mds_kwargs)
    coords = mds.fit_transform(distances.values)
    scores = None
    if ranked is not None and "anomaly_score" in ranked.columns:
        scores = ranked.reindex(X_binary.index)["anomaly_score"].astype(float).fillna(0)
    fig, ax = plt.subplots(figsize=(7, 6))
    if scores is not None:
        points = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=scores.values,
            s=70,
            cmap="viridis",
            edgecolors="black",
            linewidths=0.3,
            alpha=0.85,
        )
        fig.colorbar(points, ax=ax, label="Anomaly score")
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=70, alpha=0.85)
    labels = _ranked_label_samples(ranked, max_labels=max_labels)
    if not labels and X_binary.shape[0] <= max_labels:
        labels = set(str(s) for s in X_binary.index)
    dx, dy = _offset_for_label(coords[:, 0], coords[:, 1])
    for i, sample in enumerate(X_binary.index):
        if str(sample) in labels:
            ax.text(coords[i, 0] + dx, coords[i, 1] + dy, sample, fontsize=8)
    stress = getattr(mds, "stress_", np.nan)
    stress_label = f" (stress={stress:.3g})" if not np.isnan(stress) else ""
    ax.set_xlabel("MDS1")
    ax.set_ylabel("MDS2")
    ax.set_title(f"Jaccard MDS{stress_label}")
    ax.margins(x=0.08, y=0.08)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_anomaly_qc_plot(
    qc_df: Any,
    out_path: str,
    max_labels: int = 14,
) -> None:
    """Plot sample hash depth against cohort-unique fraction."""
    if qc_df.empty:
        return
    _, np = _load_pandas_numpy()
    plt, _ = _load_plotting()
    scores = (
        qc_df["anomaly_score"].astype(float)
        if "anomaly_score" in qc_df.columns
        else qc_df["pct_unique"].astype(float)
    )
    unique = qc_df["unique_hashes"].astype(float)
    span = float(unique.max() - unique.min())
    if span <= 1e-12:
        sizes = np.full(unique.shape, 90.0)
    else:
        sizes = 45.0 + 260.0 * ((unique - unique.min()) / span)

    fig, ax = plt.subplots(figsize=(8, 6))
    points = ax.scatter(
        qc_df["n_hashes"],
        qc_df["pct_unique"],
        c=scores,
        s=sizes,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.3,
        alpha=0.85,
    )
    fig.colorbar(
        points,
        ax=ax,
        label="Anomaly score" if "anomaly_score" in qc_df.columns else "Unique fraction",
    )
    if qc_df["n_hashes"].max() / max(1, qc_df["n_hashes"].min()) > 10:
        ax.set_xscale("log")
    labels_ordered: list[str] = []
    seen_labels: set[str] = set()
    def add_labels(samples: Iterable[Any]) -> None:
        if max_labels <= 0:
            return
        for sample in samples:
            text = str(sample)
            if text not in seen_labels and len(labels_ordered) < max_labels:
                labels_ordered.append(text)
                seen_labels.add(text)

    if "anomaly_score" in qc_df.columns:
        add_labels(
            qc_df.sort_values(
                "anomaly_score",
                ascending=False,
            ).head(max(1, max_labels // 2)).index
        )
    add_labels(
        qc_df.sort_values("n_hashes").head(max(1, max_labels // 4)).index
    )
    add_labels(
        qc_df.sort_values(
            "pct_unique",
            ascending=False,
        ).head(max(1, max_labels // 4)).index
    )
    labels = set(labels_ordered)
    dx, dy = _offset_for_label(qc_df["n_hashes"].values, qc_df["pct_unique"].values)
    for sample, row in qc_df.iterrows():
        if str(sample) in labels:
            ax.text(
                float(row["n_hashes"]) + dx,
                float(row["pct_unique"]) + dy,
                str(sample),
                fontsize=8,
            )
    ax.set_xlabel("Number of observed hashes")
    ax.set_ylabel("Fraction of hashes unique to this sample")
    ax.set_title("Anomaly QC: Hash Depth vs Cohort-Uniqueness")
    ax.grid(alpha=0.2)
    ax.margins(x=0.12, y=0.08)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_driver_hash_heatmap(
    df: Any,
    ranked: Any,
    out_path: str,
    top_hashes: int = 60,
    top_samples: int = 20,
) -> None:
    """Write a heatmap of high-difference hashes for ranked samples."""
    if df.empty or ranked.empty or top_hashes <= 0 or top_samples <= 0:
        return
    pd, np = _load_pandas_numpy()
    plt, sns = _load_plotting()
    counts = df.fillna(0).astype(float)
    ranked_samples = [s for s in ranked.index if s in counts.columns]
    if not ranked_samples:
        return
    high_samples = ranked_samples[:min(top_samples, len(ranked_samples))]
    low_samples = ranked_samples[-min(10, max(0, len(ranked_samples) - len(high_samples))):]
    sample_order = list(dict.fromkeys(high_samples + low_samples))
    foreground = counts[high_samples]
    background_cols = [c for c in counts.columns if c not in high_samples]
    if background_cols:
        background = counts[background_cols].median(axis=1)
    else:
        background = pd.Series(0.0, index=counts.index)
    driver_score = foreground.max(axis=1) - background
    selected = driver_score[driver_score > 0].sort_values(ascending=False).head(
        top_hashes
    ).index
    if len(selected) == 0:
        selected = counts.sum(axis=1).sort_values(ascending=False).head(top_hashes).index
    matrix = np.log1p(counts.loc[selected, sample_order].T)
    matrix.columns = [str(c) for c in matrix.columns]
    width = max(10.0, min(18.0, 0.18 * len(selected)))
    height = max(6.0, min(16.0, 0.28 * len(sample_order)))
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        matrix,
        cmap="magma",
        ax=ax,
        cbar_kws={"label": "log1p(count)"},
    )
    ax.set_xlabel("Driver hash")
    ax.set_ylabel("Sample")
    ax.set_title("Top Driver Hashes Across Ranked Samples")
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def compute_enrichment_plot(enrichment_df: Any, out_path: str) -> None:
    """Write a ranked -log10(p-value) plot for hash enrichment results."""
    if enrichment_df is None or enrichment_df.empty or "pvalue" not in enrichment_df:
        return
    _, np = _load_pandas_numpy()
    plt, _ = _load_plotting()
    data = enrichment_df.reset_index().copy()
    if "hash" not in data.columns:
        data = data.rename(columns={data.columns[0]: "hash"})
    data = data.sort_values("pvalue", ascending=True).reset_index(drop=True)
    pvalues = data["pvalue"].astype(float).clip(lower=np.nextafter(0, 1))
    y = -np.log10(pvalues)
    significant = (
        data["p_adj"].astype(float) <= 0.05 if "p_adj" in data.columns
        else np.zeros(len(data), dtype=bool)
    )
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = np.where(significant, "#c0392b", "#2f80ed")
    ax.scatter(range(1, len(data) + 1), y, c=colors, s=16, alpha=0.75, linewidths=0)
    ax.axhline(-math.log10(0.05), color="#52606d", linestyle="--", linewidth=1)
    if "p_adj" in data.columns and (data["p_adj"].astype(float) <= 0.05).any():
        min_sig = data.loc[data["p_adj"].astype(float) <= 0.05, "pvalue"].max()
        ax.axhline(-math.log10(max(float(min_sig), np.nextafter(0, 1))),
                   color="#c0392b", linestyle=":", linewidth=1)
    if len(data) <= 100:
        label_rows = data.head(min(8, len(data)))
        for idx, row in label_rows.iterrows():
            ax.text(
                idx + 1,
                y.iloc[idx] + 0.03,
                str(row["hash"]),
                fontsize=7,
                rotation=30,
                ha="left",
            )
    ax.set_xlabel("Hash rank by enrichment p-value")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Hash Enrichment Signal")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def combine_scores(mean_dist: Any, if_scores: Any, lof_scores: Any) -> Any:
    pd, _ = _load_pandas_numpy()
    df = pd.DataFrame({
        "mean_jaccard": mean_dist,
        "isolation_forest": if_scores,
        "lof": lof_scores,
    })
    if df.empty:
        df["combined_rank"] = []
        df["combined_score"] = []
        df["anomaly_score"] = []
        return df

    ranks = df.rank(ascending=False)
    df["combined_rank"] = ranks.mean(axis=1)
    rank_span = float(df["combined_rank"].max() - df["combined_rank"].min())
    if rank_span <= 1e-12:
        df["combined_score"] = 0.0
        df["anomaly_score"] = 0.0
    else:
        df["combined_score"] = (
            (df["combined_rank"] - df["combined_rank"].min()) / rank_span
        )
        df["anomaly_score"] = 1.0 - df["combined_score"]
    return df.sort_values("anomaly_score", ascending=False)


def top_contributing_hashes(df: Any, sample: str, top_n: int = 10) -> list[dict[str, Any]]:
    if sample not in df.columns:
        return []
    sample_counts = df[sample].fillna(0).astype(float)
    if df.shape[1] > 1:
        other = df.drop(columns=[sample]).median(axis=1).fillna(0).astype(float)
    else:
        pd, _ = _load_pandas_numpy()
        other = pd.Series(0.0, index=df.index)
    diff = sample_counts - other
    records = []
    for h in df.index:
        records.append({
            "hash": int(h),
            "count": float(sample_counts.loc[h]),
            "other_median": float(other.loc[h]),
            "diff": float(diff.loc[h]),
        })
    records.sort(key=lambda row: (row["diff"], row["count"]), reverse=True)
    return [row for row in records if row["count"] > 0 or row["diff"] > 0][:top_n]


def write_explanations(df: Any, combined: Any, out_dir: str, top_n: int = 10) -> str:
    pd, _ = _load_pandas_numpy()
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    pres = (df.fillna(0) > 0).astype(int)
    prevalence = pres.sum(axis=1)

    for sample in combined.index:
        score = (
            float(combined.loc[sample, "anomaly_score"])
            if "anomaly_score" in combined.columns else 0.0
        )
        top = top_contributing_hashes(df, sample, top_n=top_n)
        sample_presence = pres[sample] == 1
        unique_mask = sample_presence & (prevalence == 1)
        total_unique = int(unique_mask.sum())
        n_unique_top = sum(
            1 for row in top if row["other_median"] == 0 and row["count"] > 0
        )

        reasons = []
        if total_unique >= max(1, int(top_n * 0.5)):
            reasons.append(
                f"contains {total_unique} hashes not seen in other samples"
            )
        elif n_unique_top >= max(1, int(top_n * 0.2)):
            reasons.append(
                f"top-{top_n} contains {n_unique_top} unique hashes"
            )
        if score > 0.8:
            reasons.append("high anomaly score")
        if not reasons:
            reasons.append("enriched for specific hashes relative to cohort")

        rows.append({
            "sample": sample,
            "anomaly_score": score,
            "top_hashes": ";".join(str(row["hash"]) for row in top),
            "top_hash_counts": ";".join(str(int(row["count"])) for row in top),
            "top_hash_diffs": ";".join(f"{row['diff']:.1f}" for row in top),
            "total_unique_hashes": total_unique,
            "unique_in_topN": n_unique_top,
            "reason": "; ".join(reasons),
        })

    expl_df = pd.DataFrame(rows).set_index("sample")
    csv_path = os.path.join(out_dir, "explanations_by_sample.csv")
    json_path = os.path.join(out_dir, "explanations_by_sample.json")
    expl_df.to_csv(csv_path)
    with open(json_path, "w") as fh:
        json.dump(expl_df.to_dict(orient="index"), fh, indent=2)
    return csv_path


def write_html_report(
    out_dir: str,
    title: str,
    summary: dict[str, object],
    tables: dict[str, Any],
    links: dict[str, str],
) -> str:
    """Write a compact static HTML report for classification outputs."""
    os.makedirs(out_dir, exist_ok=True)
    css = """
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           margin: 2rem; color: #1f2933; }
    h1, h2 { color: #102a43; }
    table { border-collapse: collapse; width: 100%; margin: 1rem 0 2rem; }
    th, td { border: 1px solid #d9e2ec; padding: 0.45rem 0.6rem; font-size: 0.9rem; }
    th { background: #f0f4f8; text-align: left; }
    .summary { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
               gap: 0.75rem; margin: 1rem 0 2rem; }
    .metric { border: 1px solid #d9e2ec; padding: 0.75rem; border-radius: 6px; }
    .metric strong { display: block; font-size: 0.8rem; color: #52606d; }
    .plots { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
             gap: 1rem; margin: 1rem 0 2rem; }
    figure { margin: 0; border: 1px solid #d9e2ec; padding: 0.75rem; border-radius: 6px; }
    figure img { max-width: 100%; height: auto; display: block; }
    figcaption { font-size: 0.85rem; color: #52606d; margin-top: 0.5rem; }
    """
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{html.escape(title)}</title>",
        f"<style>{css}</style></head><body>",
        f"<h1>{html.escape(title)}</h1>",
        "<section class='summary'>",
    ]
    for key, value in summary.items():
        parts.append(
            "<div class='metric'>"
            f"<strong>{html.escape(str(key))}</strong>"
            f"{html.escape(str(value))}</div>"
        )
    parts.append("</section>")

    if links:
        parts.append("<h2>Artifacts</h2><ul>")
        for label, path in links.items():
            parts.append(
                f"<li><a href='{html.escape(os.path.basename(path))}'>"
                f"{html.escape(label)}</a></li>"
            )
        parts.append("</ul>")

        image_links = {
            label: path
            for label, path in links.items()
            if str(path).lower().endswith(".png")
        }
        if image_links:
            parts.append("<h2>Diagnostic Plots</h2><section class='plots'>")
            for label, path in image_links.items():
                href = html.escape(os.path.basename(path))
                parts.append(
                    "<figure>"
                    f"<a href='{href}'><img src='{href}' alt='{html.escape(label)}'></a>"
                    f"<figcaption>{html.escape(label)}</figcaption>"
                    "</figure>"
                )
            parts.append("</section>")

    for label, table in tables.items():
        parts.append(f"<h2>{html.escape(label)}</h2>")
        if hasattr(table, "to_html"):
            parts.append(table.to_html(classes="data", escape=True, border=0))
        else:
            parts.append(f"<pre>{html.escape(str(table))}</pre>")

    parts.append("</body></html>")
    report_path = os.path.join(out_dir, "outlimer_report.html")
    with open(report_path, "w") as fh:
        fh.write("\n".join(parts))
    return report_path


def write_multiqc_custom_content(
    out_dir: str,
    ranked: Any,
    score_column: str,
    filename: str = "outlimer_anomaly_scores_mqc.json",
) -> str:
    """Write MultiQC custom-content JSON for anomaly scores."""
    os.makedirs(out_dir, exist_ok=True)
    data = {}
    if score_column in ranked.columns:
        for sample, row in ranked.iterrows():
            try:
                data[str(sample)] = {"score": float(row[score_column])}
            except (TypeError, ValueError):
                continue
    payload = {
        "id": "outlimer_anomaly_scores",
        "section_name": "OutliMer anomaly scores",
        "description": "OutliMer per-sample anomaly ranking.",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "outlimer_anomaly_scores_plot",
            "title": "OutliMer anomaly scores",
            "ylab": score_column,
        },
        "data": data,
    }
    out_path = os.path.join(out_dir, filename)
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return out_path


def write_ro_crate(
    out_dir: str,
    name: str,
    inputs: dict[str, str],
    outputs: dict[str, str],
    parameters: dict[str, object],
) -> str:
    """Write a minimal RO-Crate metadata file for a run directory."""
    graph = [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
            "about": {"@id": "./"},
        },
        {
            "@id": "./",
            "@type": "Dataset",
            "name": name,
            "dateCreated": datetime.now(timezone.utc).isoformat(),
            "hasPart": [{"@id": os.path.basename(path)}
                        for path in list(inputs.values()) + list(outputs.values())
                        if path],
            "outlimerParameters": parameters,
        },
    ]
    for label, path in {**inputs, **outputs}.items():
        if not path:
            continue
        graph.append({
            "@id": os.path.basename(path),
            "@type": "File",
            "name": label,
            "contentUrl": path,
        })
    crate = {
        "@context": "https://w3id.org/ro/crate/1.1/context",
        "@graph": graph,
    }
    out_path = os.path.join(out_dir, "ro-crate-metadata.json")
    with open(out_path, "w") as fh:
        json.dump(crate, fh, indent=2)
    return out_path


def write_sparse_npz(df: Any, out_path: str) -> str:
    """Write a sparse sample-by-hash matrix as scipy .npz plus sidecar labels."""
    try:
        from scipy import sparse
    except Exception as exc:
        raise RuntimeError("Sparse NPZ export requires scipy") from exc
    matrix = sparse.csr_matrix(df.T.astype(float).values)
    sparse.save_npz(out_path, matrix)
    labels_path = out_path + ".labels.json"
    with open(labels_path, "w") as fh:
        json.dump({
            "samples": [str(s) for s in df.columns],
            "hashes": [str(h) for h in df.index],
            "orientation": "samples_by_hashes",
        }, fh, indent=2)
    return out_path


def compare_union_matrices(baseline_csv: str, query_csv: str) -> Any:
    """Compare two union matrices sample-by-sample against a baseline hash set."""
    pd, _ = _load_pandas_numpy()
    baseline = load_union_csv(baseline_csv)
    query = load_union_csv(query_csv)
    baseline_hashes = set(int(h) for h in baseline.index[(baseline > 0).any(axis=1)])
    rows = []
    for sample in query.columns:
        sample_hashes = set(int(h) for h in query.index[query[sample] > 0])
        n_total = len(sample_hashes)
        n_in = len(sample_hashes & baseline_hashes)
        pct = n_in / n_total if n_total else 0.0
        rows.append({
            "sample": sample,
            "n_hashes": n_total,
            "n_in_baseline": n_in,
            "pct_in_baseline": pct,
            "n_new_hashes": n_total - n_in,
        })
    return pd.DataFrame(rows).set_index("sample").sort_values(
        ["pct_in_baseline", "n_new_hashes"],
        ascending=[True, False],
    )


def _validate_sample_subset(
    label: str,
    requested: Sequence[str],
    available: set[str],
) -> list[str]:
    unknown = sorted(set(requested) - available)
    if unknown:
        raise ValueError(f"Unknown {label} sample(s): {', '.join(unknown)}")
    return list(dict.fromkeys(requested))


def _benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda idx: pvalues[idx])
    adjusted = [1.0] * n
    running = 1.0
    for sorted_pos in range(n - 1, -1, -1):
        idx = order[sorted_pos]
        rank = sorted_pos + 1
        running = min(running, float(pvalues[idx]) * n / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def compute_enrichment(
    sample_hash_counts: dict[str, dict[int, int]],
    foreground_samples: Sequence[str],
    background_samples: Sequence[str] | None = None,
) -> Any:
    """Compute foreground-vs-background per-hash Fisher exact enrichment."""
    pd, _ = _load_pandas_numpy()
    fisher_exact = _load_fisher_exact()

    available = set(sample_hash_counts)
    foreground = _validate_sample_subset("foreground", foreground_samples, available)
    if background_samples is None:
        background = sorted(available - set(foreground))
    else:
        background = _validate_sample_subset("background", background_samples, available)

    overlap = set(foreground) & set(background)
    if overlap:
        raise ValueError(
            "Foreground and background samples must be disjoint: "
            + ", ".join(sorted(overlap))
        )
    if not foreground:
        raise ValueError("At least one foreground sample is required")
    if not background:
        raise ValueError("At least one background sample is required")

    all_hashes = sorted({
        int(h)
        for counts in sample_hash_counts.values()
        for h, c in counts.items()
        if int(c) > 0
    })
    rows = []
    for h in all_hashes:
        fg_present = sum(1 for s in foreground if sample_hash_counts[s].get(h, 0) > 0)
        bg_present = sum(1 for s in background if sample_hash_counts[s].get(h, 0) > 0)
        table = [
            [fg_present, len(foreground) - fg_present],
            [bg_present, len(background) - bg_present],
        ]
        oddsratio, pvalue = fisher_exact(table)
        rows.append({
            "hash": h,
            "pvalue": float(pvalue),
            "oddsratio": float(oddsratio),
            "count_in_foreground": fg_present,
            "count_in_background": bg_present,
            "n_foreground": len(foreground),
            "n_background": len(background),
            "total_count": sum(
                int(counts.get(h, 0)) for counts in sample_hash_counts.values()
            ),
        })

    res = pd.DataFrame(rows)
    if res.empty:
        return pd.DataFrame(columns=[
            "pvalue", "oddsratio", "count_in_foreground",
            "count_in_background", "n_foreground", "n_background",
            "total_count", "p_adj",
        ])
    res["p_adj"] = _benjamini_hochberg(res["pvalue"].fillna(1.0).tolist())
    return res.set_index("hash").sort_values(["p_adj", "pvalue"])


def _counts_from_union(df: Any) -> dict[str, dict[int, int]]:
    return {
        sample: {
            int(h): int(v)
            for h, v in df[sample].fillna(0).astype(int).to_dict().items()
            if int(v) > 0
        }
        for sample in df.columns
    }


def _default_report_path(out_dir: str, report_csv: str | None) -> str:
    if report_csv:
        return report_csv
    candidate = os.path.join(out_dir, "kmer_report.csv")
    if os.path.exists(candidate):
        return candidate
    return "kmer_report.csv"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Classify OutliMer k-mer hash matrices and reports."
    )
    parser.add_argument(
        "--union-csv",
        help="Top-union summary CSV (rows=hash, columns=samples)",
    )
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "--top-features",
        type=int,
        default=2000,
        help="Max number of hash features to keep (0 = no limit)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum sample prevalence for a hash to be retained",
    )
    parser.add_argument(
        "--contamination",
        type=float,
        default=0.05,
        help="Expected anomaly fraction for IF/LOF; must be >0 and <=0.5",
    )
    parser.add_argument("--no-plots", action="store_true", help="Do not write plots")
    parser.add_argument(
        "--plot-labels",
        type=int,
        default=12,
        help="Maximum high-priority sample labels on ordination/QC plots",
    )
    parser.add_argument(
        "--plot-top-hashes",
        type=int,
        default=60,
        help="Number of hashes to show in the driver-hash heatmap",
    )
    parser.add_argument(
        "--plot-top-samples",
        type=int,
        default=20,
        help="Number of top-ranked samples to show in the driver-hash heatmap",
    )
    parser.add_argument(
        "--report-csv",
        default=None,
        help="OutliMer kmer_report.csv for report mode",
    )
    parser.add_argument(
        "--mode",
        choices=["union", "report", "compare", "auto"],
        default="auto",
        help="Use union features, report heuristics, compare runs, or auto-detect inputs",
    )
    parser.add_argument(
        "--explain-top-n",
        type=int,
        default=10,
        help="Top-N hashes to include in explanations",
    )
    parser.add_argument(
        "--explain-output",
        default=None,
        help="Directory for explanations (defaults to --out-dir)",
    )
    parser.add_argument(
        "--foreground-samples",
        help="Comma-separated foreground sample names for enrichment",
    )
    parser.add_argument(
        "--background-samples",
        help="Comma-separated background sample names for enrichment",
    )
    parser.add_argument(
        "--foreground-query",
        help="Pandas metadata query selecting foreground samples",
    )
    parser.add_argument(
        "--background-query",
        help="Pandas metadata query selecting background samples",
    )
    parser.add_argument(
        "--baseline-union",
        help="Baseline union CSV for --mode compare",
    )
    parser.add_argument(
        "--query-union",
        help="Query union CSV for --mode compare",
    )
    parser.add_argument(
        "--sparse-npz",
        nargs="?",
        const="feature_matrix.npz",
        help="Write sparse sample-by-hash matrix as .npz (optional filename)",
    )
    parser.add_argument(
        "--multiqc",
        action="store_true",
        help="Write MultiQC custom-content JSON",
    )
    parser.add_argument(
        "--ro-crate",
        action="store_true",
        help="Write minimal ro-crate-metadata.json provenance",
    )
    parser.add_argument("--metadata", help="Sample metadata CSV/TSV")
    parser.add_argument(
        "--metadata-sample-column",
        default="sample",
        help="Column in --metadata containing sample names (default: sample)",
    )
    parser.add_argument(
        "--no-html-report",
        action="store_true",
        help="Do not write outlimer_report.html",
    )
    return parser


def _select_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    if args.mode == "compare":
        if not args.baseline_union or not args.query_union:
            parser.error("--mode compare requires --baseline-union and --query-union")
        return "compare"
    if args.mode == "union":
        if not args.union_csv:
            parser.error("--mode union requires --union-csv")
        return "union"
    if args.mode == "report":
        return "report"
    if args.union_csv and os.path.exists(args.union_csv):
        return "union"
    return "report"


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        _validate_contamination(args.contamination)
    except ValueError as exc:
        parser.error(str(exc))
    if args.top_features < 0:
        parser.error("--top-features must be >= 0")
    if args.min_samples < 1:
        parser.error("--min-samples must be >= 1")
    if args.explain_top_n < 1:
        parser.error("--explain-top-n must be >= 1")
    if args.plot_labels < 0:
        parser.error("--plot-labels must be >= 0")
    if args.plot_top_hashes < 0:
        parser.error("--plot-top-hashes must be >= 0")
    if args.plot_top_samples < 0:
        parser.error("--plot-top-samples must be >= 0")
    if (args.foreground_query or args.background_query) and not args.metadata:
        parser.error("--foreground-query/--background-query require --metadata")

    os.makedirs(args.out_dir, exist_ok=True)
    mode_used = _select_mode(args, parser)
    metadata = None
    if args.metadata:
        try:
            metadata = load_sample_metadata(args.metadata, args.metadata_sample_column)
        except Exception as exc:
            parser.error(f"Failed to load --metadata {args.metadata}: {exc}")

    if mode_used == "compare":
        try:
            compared = compare_union_matrices(args.baseline_union, args.query_union)
            compared = attach_metadata(compared, metadata)
        except Exception as exc:
            print(f"Run comparison failed: {exc}", file=sys.stderr)
            return 1
        compare_path = os.path.join(args.out_dir, "run_comparison.csv")
        compared.to_csv(compare_path)
        print(f"Wrote run comparison to {compare_path}")
        links = {"Run comparison CSV": compare_path}
        if args.multiqc:
            mqc_scores = compared.copy()
            mqc_scores["anomaly_score"] = 1.0 - mqc_scores["pct_in_baseline"]
            mqc_path = write_multiqc_custom_content(
                args.out_dir,
                mqc_scores,
                "anomaly_score",
                filename="outlimer_run_comparison_mqc.json",
            )
            links["MultiQC custom content"] = mqc_path
            print(f"Wrote MultiQC custom content to {mqc_path}")
        if not args.no_html_report:
            html_path = write_html_report(
                args.out_dir,
                "OutLiMer Run Comparison",
                {
                    "mode": "compare",
                    "query samples": compared.shape[0],
                    "baseline": args.baseline_union,
                    "query": args.query_union,
                },
                {"Most novel query samples": compared.head(25)},
                links,
            )
            print(f"Wrote HTML report to {html_path}")
        if args.ro_crate:
            crate_path = write_ro_crate(
                args.out_dir,
                "OutliMer run comparison",
                {"baseline union": args.baseline_union, "query union": args.query_union},
                links,
                {"mode": "compare"},
            )
            print(f"Wrote RO-Crate metadata to {crate_path}")

    elif mode_used == "union":
        plot_links: dict[str, str] = {}
        try:
            df = load_union_csv(args.union_csv)
            X, X_binary, X_log = prepare_feature_matrix(
                df,
                top_M=args.top_features,
                min_samples=args.min_samples,
            )
            mean_dist = compute_mean_jaccard_distance(X_binary)
            if_scores = compute_isolation_forest(
                X_log,
                contamination=args.contamination,
            )
            lof_scores = compute_lof(X_log, contamination=args.contamination)
            combined = combine_scores(mean_dist, if_scores, lof_scores)
            if not args.no_plots:
                qc_metrics = compute_sample_qc_metrics(df, combined)
                pca_path = os.path.join(args.out_dir, "pca_samples.png")
                compute_pca_plot(
                    X_log,
                    pca_path,
                    ranked=combined,
                    max_labels=args.plot_labels,
                )
                plot_links["PCA plot"] = pca_path

                mds_path = os.path.join(args.out_dir, "jaccard_mds.png")
                compute_jaccard_mds_plot(
                    X_binary,
                    mds_path,
                    ranked=combined,
                    max_labels=args.plot_labels,
                )
                plot_links["Jaccard MDS plot"] = mds_path

                qc_path = os.path.join(args.out_dir, "anomaly_qc_scatter.png")
                compute_anomaly_qc_plot(
                    qc_metrics,
                    qc_path,
                    max_labels=args.plot_labels,
                )
                plot_links["Anomaly QC scatter"] = qc_path

                dendrogram_path = os.path.join(args.out_dir, "dendrogram.png")
                compute_dendrogram(X_binary, dendrogram_path)
                plot_links["Jaccard dendrogram"] = dendrogram_path

                cluster_path = os.path.join(args.out_dir, "jaccard_cluster_heatmap.png")
                compute_clustered_jaccard_heatmap(X_binary, cluster_path)
                plot_links["Clustered Jaccard heatmap"] = cluster_path

                driver_path = os.path.join(args.out_dir, "driver_hash_heatmap.png")
                compute_driver_hash_heatmap(
                    df,
                    combined,
                    driver_path,
                    top_hashes=args.plot_top_hashes,
                    top_samples=args.plot_top_samples,
                )
                plot_links["Driver hash heatmap"] = driver_path
            combined = attach_metadata(combined, metadata)
        except Exception as exc:
            print(f"Union-mode classification failed: {exc}", file=sys.stderr)
            return 1

        links = dict(plot_links)
        if args.sparse_npz:
            sparse_name = args.sparse_npz
            sparse_path = sparse_name if os.path.isabs(sparse_name) else os.path.join(
                args.out_dir, sparse_name)
            try:
                write_sparse_npz(df, sparse_path)
                links["Sparse feature matrix"] = sparse_path
                print(f"Wrote sparse feature matrix to {sparse_path}")
            except Exception as exc:
                print(f"Failed sparse matrix export: {exc}", file=sys.stderr)

        out_report = os.path.join(args.out_dir, "outliers_report.csv")
        combined.to_csv(out_report, index=True)
        print(f"Wrote outlier report to {out_report}")
        links["Outlier report CSV"] = out_report

        topn = max(1, int(math.ceil(0.05 * combined.shape[0]))) if not combined.empty else 0
        topn_path = os.path.join(args.out_dir, "top_anomalies.csv")
        combined.head(topn).to_csv(topn_path)
        print(f"Wrote top {topn} anomalies to {topn_path}")
        links["Top anomalies CSV"] = topn_path

        explain_outdir = args.explain_output or args.out_dir
        try:
            expl_path = write_explanations(
                df,
                combined,
                explain_outdir,
                top_n=args.explain_top_n,
            )
            print(f"Wrote per-sample explanations to {expl_path}")
            links["Sample explanations"] = expl_path
        except Exception as exc:
            print(f"Failed to write explanations: {exc}", file=sys.stderr)

        try:
            sample_hash_counts = _counts_from_union(df)
            available = set(sample_hash_counts)
            foreground = _parse_sample_list(args.foreground_samples)
            background = _parse_sample_list(args.background_samples)
            if args.foreground_query:
                foreground = select_samples_by_query(metadata, args.foreground_query)
            if args.background_query:
                background = select_samples_by_query(metadata, args.background_query)
            if foreground is None:
                foreground = combined.head(topn).index.tolist()
            if background is None:
                background = sorted(available - set(foreground))
            if foreground and background:
                enr = compute_enrichment(sample_hash_counts, foreground, background)
                enr_path = os.path.join(explain_outdir, "hash_enrichment.csv")
                enr.to_csv(enr_path)
                links["Hash enrichment CSV"] = enr_path
                print(f"Wrote enrichment results to {enr_path}")
                if not args.no_plots:
                    enr_plot = os.path.join(explain_outdir, "hash_enrichment_plot.png")
                    compute_enrichment_plot(enr, enr_plot)
                    links["Hash enrichment plot"] = enr_plot
            else:
                print(
                    "Skipping enrichment: need at least one foreground and "
                    "one background sample",
                    file=sys.stderr,
                )
        except Exception as exc:
            print(f"Failed enrichment computation: {exc}", file=sys.stderr)

        if args.multiqc:
            mqc_path = write_multiqc_custom_content(
                args.out_dir, combined, "anomaly_score")
            links["MultiQC custom content"] = mqc_path
            print(f"Wrote MultiQC custom content to {mqc_path}")

        if not args.no_html_report:
            report_path = write_html_report(
                args.out_dir,
                "OutliMer Classification Report",
                {
                    "mode": "union",
                    "samples": combined.shape[0],
                    "features": df.shape[0],
                    "top anomalies": topn,
                },
                {
                    "Top anomalies": combined.head(max(topn, 1)),
                },
                links,
            )
            print(f"Wrote HTML report to {report_path}")
            links["HTML report"] = report_path
        if args.ro_crate:
            crate_path = write_ro_crate(
                args.out_dir,
                "OutliMer classification",
                {"union csv": args.union_csv, "metadata": args.metadata or ""},
                links,
                {
                    "mode": "union",
                    "top_features": args.top_features,
                    "min_samples": args.min_samples,
                    "contamination": args.contamination,
                },
            )
            print(f"Wrote RO-Crate metadata to {crate_path}")

    else:
        report_path = _default_report_path(args.out_dir, args.report_csv)
        try:
            report_df = load_report_csv(report_path)
            explained = explain_report_based(report_df)
            explained = attach_metadata(explained, metadata)
        except Exception as exc:
            print(f"Report-mode classification failed: {exc}", file=sys.stderr)
            return 1

        links = {}
        out_report = os.path.join(args.out_dir, "outliers_report_from_reportcsv.csv")
        explained.to_csv(out_report)
        print(f"Wrote outlier report (from OutliMer report) to {out_report}")
        links["Report-mode CSV"] = out_report
        text_path = os.path.join(args.out_dir, "outlier_explanations.txt")
        with open(text_path, "w") as tf:
            for sample in explained.index:
                row = explained.loc[sample]
                tf.write(
                    f"{sample}: score={row['reason_score']:.4f}; "
                    f"reasons={row['reasons']}; pct_in_db={row['pct_in_db']:.3f}; "
                    f"n_hashes={int(row['n_hashes'])}; "
                    f"n_new_hashes={int(row['n_new_hashes'])}\n"
                )
        print(f"Wrote textual explanations to {text_path}")
        links["Text explanations"] = text_path

        if args.multiqc:
            mqc_path = write_multiqc_custom_content(
                args.out_dir,
                explained.rename(columns={"reason_score": "anomaly_score"}),
                "anomaly_score",
            )
            links["MultiQC custom content"] = mqc_path
            print(f"Wrote MultiQC custom content to {mqc_path}")

        if not args.no_html_report:
            html_path = write_html_report(
                args.out_dir,
                "OutliMer Report-Mode Summary",
                {
                    "mode": "report",
                    "samples": explained.shape[0],
                    "source report": report_path,
                },
                {
                    "Ranked samples": explained.head(25),
                },
                links,
            )
            print(f"Wrote HTML report to {html_path}")
            links["HTML report"] = html_path
        if args.ro_crate:
            crate_path = write_ro_crate(
                args.out_dir,
                "OutliMer report-mode classification",
                {"report csv": report_path, "metadata": args.metadata or ""},
                links,
                {"mode": "report"},
            )
            print(f"Wrote RO-Crate metadata to {crate_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
