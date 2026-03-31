import argparse
import os
import sys
import math
from typing import Tuple
from typing import List, Dict

try:
	import pandas as pd
	import numpy as np
	from sklearn.decomposition import PCA
	from sklearn.ensemble import IsolationForest
	from sklearn.neighbors import LocalOutlierFactor
	from scipy.spatial.distance import pdist, squareform
	from scipy.cluster.hierarchy import linkage, dendrogram
	from scipy.stats import fisher_exact, chi2_contingency
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import seaborn as sns
except Exception as e:
	print("Missing required Python packages. Please install requirements.txt in this directory.", file=sys.stderr)
	raise


def load_union_csv(path: str) -> pd.DataFrame:
	if not os.path.exists(path):
		raise FileNotFoundError(path)
	df = pd.read_csv(path, index_col=0)
	# drop 'sequence' and 'total_count' if present
	for c in ('sequence', 'total_count'):
		if c in df.columns:
			df = df.drop(columns=[c])
	# ensure numeric
	df = df.fillna(0)
	# columns are samples? top_union assumed: rows=hash, cols = total_count, samples...
	# after drop, df columns should be samples
	return df


def load_report_csv(path: str) -> pd.DataFrame:
	"""Load OutliMer report.csv which contains per-sample summary fields.

	Expected columns: sample, n_hashes, n_in_db, pct_in_db, n_new_hashes
	Returns DataFrame indexed by sample.
	"""
	if not os.path.exists(path):
		raise FileNotFoundError(path)
	df = pd.read_csv(path, index_col=0)
	# normalize column names if needed
	cols = [c.lower() for c in df.columns]
	df.columns = cols
	# ensure expected columns exist
	expected = {'n_hashes', 'n_in_db', 'pct_in_db', 'n_new_hashes'}
	if not expected.issubset(set(df.columns)):
		# try to detect variants
		raise ValueError(f"Report CSV missing expected columns: {expected} (found {set(df.columns)})")
	# coerce types
	df['n_hashes'] = df['n_hashes'].astype(float)
	df['n_in_db'] = df['n_in_db'].astype(float)
	df['pct_in_db'] = df['pct_in_db'].astype(float)
	df['n_new_hashes'] = df['n_new_hashes'].astype(float)
	return df


def explain_report_based(report_df: pd.DataFrame) -> pd.DataFrame:
	"""Generate simple heuristic explanations for outliers based on OutliMer report.

	Returns DataFrame with columns: sample, reason_score, reasons_text, pct_in_db, n_hashes, n_new_hashes
	Higher reason_score -> more anomalous.
	"""
	df = report_df.copy()
	# compute normalized components
	# pct component: low pct_in_db => anomalous
	pct_comp = 1.0 - df['pct_in_db']
	# new hashes component
	max_new = max(1.0, df['n_new_hashes'].max())
	new_comp = df['n_new_hashes'] / max_new
	# size penalty: very small samples can be noisy
	max_size = max(1.0, df['n_hashes'].max())
	size_comp = 1.0 - (df['n_hashes'] / max_size)

	# weights
	w_pct, w_new, w_size = 0.6, 0.3, 0.1
	score = (w_pct * pct_comp) + (w_new * new_comp) + (w_size * size_comp)
	# normalize to 0-1
	score = (score - score.min()) / (score.max() - score.min() + 1e-12)

	reasons = []
	for s, row in df.iterrows():
		r = []
		if row['pct_in_db'] < 0.2:
			r.append(f"low overlap with DB (pct_in_db={row['pct_in_db']:.3f})")
		elif row['pct_in_db'] < 0.5:
			r.append(f"moderate overlap with DB (pct_in_db={row['pct_in_db']:.3f})")
		else:
			r.append(f"high overlap with DB (pct_in_db={row['pct_in_db']:.3f})")

		if row['n_new_hashes'] > (df['n_new_hashes'].median() * 2):
			r.append(f"many novel hashes (n_new_hashes={int(row['n_new_hashes'])})")
		elif row['n_new_hashes'] > 0:
			r.append(f"some novel hashes (n_new_hashes={int(row['n_new_hashes'])})")

		if row['n_hashes'] < max(50, int(max_size * 0.01)):
			r.append(f"very small sample (n_hashes={int(row['n_hashes'])}); may be noisy")

		reasons.append('; '.join(r))

	out = pd.DataFrame({
		'sample': df.index,
		'pct_in_db': df['pct_in_db'].values,
		'n_hashes': df['n_hashes'].values,
		'n_new_hashes': df['n_new_hashes'].values,
		'reason_score': score.values,
		'reasons': reasons
	}).set_index('sample')
	out = out.sort_values('reason_score', ascending=False)
	return out


def prepare_feature_matrix(df: pd.DataFrame, top_M:int=2000, min_samples:int=1) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
	# transpose to samples x features
	X = df.T.astype(float).fillna(0)
	# filter features by prevalence
	prevalence = X.gt(0).sum(axis=0)
	keep_mask = prevalence >= min_samples
	X = X.loc[:, keep_mask]
	# select top M features by total count
	if top_M and X.shape[1] > top_M:
		totals = X.sum(axis=0)
		top_feats = totals.sort_values(ascending=False).head(top_M).index
		X = X.loc[:, top_feats]
	X_binary = (X > 0).astype(int)
	X_log = np.log1p(X)
	return X, X_binary, X_log


def compute_pca_plot(X_log: pd.DataFrame, out_path: str) -> None:
	n_samples, n_features = X_log.shape
	if n_samples < 2 or n_features < 2:
		print(f"Skipping PCA plot: need >=2 samples and >=2 features (have {n_samples} samples, {n_features} features)")
		return
	n_comp = min(4, n_samples, n_features)
	if n_comp < 2:
		print(f"Skipping PCA plot: not enough components (n_comp={n_comp})")
		return
	pca = PCA(n_components=n_comp)
	Z = pca.fit_transform(X_log.values)
	fig, ax = plt.subplots(figsize=(7,6))
	sns.scatterplot(x=Z[:,0], y=Z[:,1], s=60, ax=ax)
	for i, sname in enumerate(X_log.index):
		ax.text(Z[i,0], Z[i,1], sname, fontsize=8)
	ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
	ax.set_title('PCA (log1p counts)')
	fig.tight_layout()
	fig.savefig(out_path, dpi=150)
	plt.close(fig)


def compute_mean_jaccard_distance(X_binary: pd.DataFrame) -> pd.Series:
	# compute pairwise Jaccard distances between samples (rows)
	arr = X_binary.values
	# if only one sample, return zero
	if arr.shape[0] <= 1:
		return pd.Series([0.0], index=X_binary.index)
	dist_vec = pdist(arr, metric='jaccard')
	dist_mat = squareform(dist_vec)
	mean_dist = dist_mat.mean(axis=1)
	return pd.Series(mean_dist, index=X_binary.index)


def compute_isolation_forest(X_log: pd.DataFrame, contamination: float=0.05) -> pd.Series:
	if X_log.shape[0] <= 1:
		return pd.Series([0.0], index=X_log.index)
	clf = IsolationForest(n_estimators=200, contamination=contamination, random_state=0)
	clf.fit(X_log.values)
	# score_samples -> higher is better (inlier), we invert to make higher = more anomalous
	scores = -clf.score_samples(X_log.values)
	return pd.Series(scores, index=X_log.index)


def compute_lof(X_log: pd.DataFrame, contamination: float=0.05) -> pd.Series:
	if X_log.shape[0] <= 1:
		return pd.Series([0.0], index=X_log.index)
	n_neighbors = min(20, max(2, X_log.shape[0]-1))
	lof = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=contamination, metric='euclidean')
	y_pred = lof.fit_predict(X_log.values)
	scores = -lof.negative_outlier_factor_
	return pd.Series(scores, index=X_log.index)


def compute_dendrogram(X_binary: pd.DataFrame, out_path: str) -> None:
	# create dendrogram from Jaccard distances
	arr = X_binary.values
	if arr.shape[0] <= 1:
		return
	dist_vec = pdist(arr, metric='jaccard')
	Z = linkage(dist_vec, method='average')
	fig = plt.figure(figsize=(10, 6))
	dendrogram(Z, labels=X_binary.index.tolist(), leaf_rotation=90)
	plt.title('Hierarchical clustering (Jaccard)')
	plt.tight_layout()
	fig.savefig(out_path, dpi=150)
	plt.close(fig)


def combine_scores(mean_dist: pd.Series, if_scores: pd.Series, lof_scores: pd.Series) -> pd.DataFrame:
	df = pd.DataFrame({'mean_jaccard': mean_dist, 'isolation_forest': if_scores, 'lof': lof_scores})
	# rank each (higher score -> rank 1 is most anomalous, so use rank(ascending=False))
	ranks = df.rank(ascending=False)
	df['combined_rank'] = ranks.mean(axis=1)
	# normalised combined score
	df['combined_score'] = (df['combined_rank'] - df['combined_rank'].min()) / (df['combined_rank'].max() - df['combined_rank'].min() + 1e-12)
	# lower combined_rank means more anomalous (we used descending rank), invert to make higher=more anomalous
	df['anomaly_score'] = 1.0 - df['combined_score']
	df = df.sort_values('anomaly_score', ascending=False)
	return df


def top_contributing_hashes(df: pd.DataFrame, sample: str, top_n: int = 10) -> List[Dict]:
	"""Return top contributing hashes for a sample.

	Strategy:
	- compute other_median (median across other samples)
	- compute diff = sample_count - other_median
	- prefer hashes with large positive diff; also mark unique hashes where other_median==0
	Returns list of dicts: {'hash': hash, 'count': c, 'other_median': m, 'diff': d}
	"""
	if sample not in df.columns:
		return []
	sample_counts = df[sample].fillna(0).astype(float)
	if df.shape[1] > 1:
		other = df.drop(columns=[sample]).median(axis=1).fillna(0).astype(float)
	else:
		other = pd.Series(0.0, index=df.index)

	diff = sample_counts - other
	recs = []
	for h in df.index:
		c = float(sample_counts.loc[h])
		m = float(other.loc[h]) if h in other.index else 0.0
		d = float(diff.loc[h])
		recs.append({'hash': int(h), 'count': c, 'other_median': m, 'diff': d})

	# sort primarily by diff desc, then by count desc
	recs = sorted(recs, key=lambda r: (r['diff'], r['count']), reverse=True)
	# filter out zeros
	recs = [r for r in recs if (r['count'] > 0 or r['diff'] > 0)]
	return recs[:top_n]


def write_explanations(df: pd.DataFrame, combined: pd.DataFrame, out_dir: str, top_n: int = 10) -> str:
	"""Write per-sample explanations CSV and JSON. Returns base path written (CSV path)."""
	os.makedirs(out_dir, exist_ok=True)
	rows = []
	# Precompute presence/absence matrix for full df (rows=hash, cols=samples)
	pres = (df.fillna(0) > 0).astype(int)
	prevalence = pres.sum(axis=1)  # number of samples containing each hash

	for sample in combined.index:
		score = float(combined.loc[sample, 'anomaly_score']) if 'anomaly_score' in combined.columns else float(0.0)
		top = top_contributing_hashes(df, sample, top_n=top_n)
		top_hashes = ';'.join([str(r['hash']) for r in top])
		top_counts = ';'.join([str(int(r['count'])) for r in top])
		top_diffs = ';'.join([f"{r['diff']:.1f}" for r in top])

		# total unique hashes for this sample across the entire feature set:
		# unique means hash present in this sample and absent from all other samples
		sample_presence = pres[sample] == 1
		unique_mask = sample_presence & (prevalence == 1)
		total_unique = int(unique_mask.sum())

		# unique among top-N reported features
		n_unique_top = sum(1 for r in top if r['other_median'] == 0 and r['count'] > 0)

		reasons = []
		# Reasoning: prefer to report total unique count (full set), not just top-N
		if total_unique >= max(1, int(len(df.columns) * 0.01)) or total_unique >= max(1, int(top_n * 0.5)):
			reasons.append(f"contains {total_unique} hashes not seen in other samples (total unique)")
		elif n_unique_top >= max(1, int(top_n * 0.2)):
			reasons.append(f"top-{top_n} contains {n_unique_top} unique hashes relative to cohort")

		if score > 0.8:
			reasons.append("high anomaly score")
		if not reasons:
			reasons.append("enriched for specific hashes relative to cohort")

		rows.append({
			'sample': sample,
			'anomaly_score': score,
			'top_hashes': top_hashes,
			'top_hash_counts': top_counts,
			'top_hash_diffs': top_diffs,
			'total_unique_hashes': total_unique,
			'unique_in_topN': n_unique_top,
			'reason': '; '.join(reasons),
		})

	expl_df = pd.DataFrame(rows).set_index('sample')
	csv_path = os.path.join(out_dir, 'explanations_by_sample.csv')
	json_path = os.path.join(out_dir, 'explanations_by_sample.json')
	expl_df.to_csv(csv_path)
	return csv_path


def compute_enrichment(sample_hash_counts: Dict[str, Dict[int, int]], background_samples: List[str] | None = None) -> pd.DataFrame:
	"""Compute per-hash enrichment across samples using Fisher's exact test.

	Returns DataFrame with columns: hash, pvalue, oddsratio, count_in_samples, count_in_background, total_in_samples, total_in_background
	"""
	# Build total counts per hash across all samples
	total_counts: Dict[int, int] = {}
	sample_presence: Dict[int, int] = {}
	samples = sorted(sample_hash_counts.keys())
	for s in samples:
		counts = sample_hash_counts[s]
		for h, c in counts.items():
			ih = int(h)
			total_counts[ih] = total_counts.get(ih, 0) + int(c)
			sample_presence[ih] = sample_presence.get(ih, 0) + (1 if int(c) > 0 else 0)

	# choose background set (if None, use all samples)
	if background_samples is None:
		background_samples = samples

	# precompute totals
	samples_set = set(samples)
	bg_set = set(background_samples)
	n_samples = len(samples_set)
	n_bg = len(bg_set)

	# counts per hash in sample set and background set
	rows = []
	for h, tot in sorted(total_counts.items(), key=lambda x: -x[1]):
		# count number of samples with this hash in sample set and background set
		count_in_samples = sum(1 for s in samples if int(h) in sample_hash_counts[s] and sample_hash_counts[s].get(int(h), 0) > 0)
		count_in_bg = sum(1 for s in background_samples if int(h) in sample_hash_counts[s] and sample_hash_counts[s].get(int(h), 0) > 0)

		# contingency table: [[in_samples, not_in_samples], [in_bg, not_in_bg]]
		a = count_in_samples
		b = n_samples - a
		c = count_in_bg
		d = n_bg - c
		# Fisher's exact (handle degenerate cases)
		try:
			oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
		except Exception:
			oddsratio, pvalue = float('nan'), 1.0

		rows.append({'hash': int(h), 'pvalue': pvalue, 'oddsratio': oddsratio, 'count_in_samples': a, 'count_in_background': c, 'total_count': int(tot)})

	res = pd.DataFrame(rows).set_index('hash')
	# multiple testing correction (Benjamini-Hochberg)
	from math import isnan
	pvals = res['pvalue'].fillna(1.0).values
	n = len(pvals)
	order = sorted(range(n), key=lambda i: pvals[i])
	bh = [0.0] * n
	for i, idx in enumerate(order):
		bh[idx] = pvals[idx] * n / (i+1)
	# ensure monotonic
	for i in range(n-2, -1, -1):
		if bh[i] > bh[i+1]:
			bh[i] = bh[i+1]
	res['p_adj'] = [min(1.0, x) for x in bh]
	res = res.sort_values('p_adj')
	return res





def main(argv=None):
	ap = argparse.ArgumentParser(description='Simple outlier detection on k-mer union summary')
	ap.add_argument('--union-csv', help='Top-union summary CSV (rows=hash, cols=samples)',
					required=True)
	ap.add_argument('--out-dir', help='Output directory',
					required=True)
	ap.add_argument('--top-features', type=int, default=2000, help='Max number of features (hashes) to keep')
	ap.add_argument('--min-samples', type=int, default=1, help='Minimum sample prevalence for a hash to be retained')
	ap.add_argument('--contamination', type=float, default=0.05, help='Expected fraction of anomalies (for IF/LOF)')
	ap.add_argument('--no-plots', action='store_true', help='Do not write plots')
	ap.add_argument('--report-csv', default='kmer_report.csv', help='OutliMer report CSV fallback (used when --union-csv missing)')
	ap.add_argument('--mode', choices=['union', 'report', 'auto'], default='auto', help='Mode: use union features (union), use OutliMer report (report), or auto to pick available')
	ap.add_argument('--explain-top-n', type=int, default=10, help='Top-N hashes to include in explanations')
	ap.add_argument('--explain-output', default=None, help='Directory to write explanations (defaults to out-dir)')
	args = ap.parse_args(argv)

	os.makedirs(args.out_dir, exist_ok=True)
	df = None
	mode_used = None
	seq_map = None
	report_df = None
	if args.mode in ('union', 'auto'):
		try:
			# load union summary; primary df used for features will drop 'sequence' if present
			df = load_union_csv(args.union_csv)
			# attempt to load sequence map (if union CSV included sequences)
			seq_map = None
			try:
				raw_df = __import__('pandas').read_csv(args.union_csv, index_col=0)
				if 'sequence' in raw_df.columns:
					# build mapping from hash (index) to sequence (first sequence if multiple)
					seq_map = {}
					for idx, row in raw_df.iterrows():
						try:
							h = int(idx)
							seq_map[h] = str(row['sequence'])
						except Exception:
							# skip non-int indices
							continue
			except Exception:
				seq_map = None
			mode_used = 'union'
		except Exception:
			df = None
			if args.mode == 'union':
				raise
	if df is None and args.mode in ('report', 'auto'):
		try:
			report_df = load_report_csv(args.report_csv)
			mode_used = 'report'
		except Exception as e:
			if args.mode == 'report':
				raise
			raise RuntimeError(f"Failed to load union or report CSV: {e}")

	if mode_used == 'union':
		X, X_binary, X_log = prepare_feature_matrix(df, top_M=args.top_features, min_samples=args.min_samples)

		# PCA plot
		pca_path = os.path.join(args.out_dir, 'pca_samples.png')
		dend_path = os.path.join(args.out_dir, 'dendrogram.png')
		if not args.no_plots:
			compute_pca_plot(X_log, pca_path)
			compute_dendrogram(X_binary, dend_path)

		# compute distances and anomaly scores regardless of plotting preference
		mean_dist = compute_mean_jaccard_distance(X_binary)
		if_scores = compute_isolation_forest(X_log, contamination=args.contamination)
		lof_scores = compute_lof(X_log, contamination=args.contamination)

		combined = combine_scores(mean_dist, if_scores, lof_scores)

		# write report
		out_report = os.path.join(args.out_dir, 'outliers_report.csv')
		combined.to_csv(out_report, index=True)
		print(f'Wrote outlier report to {out_report}')
		if not args.no_plots:
			print(f'Wrote plots: {pca_path}, {dend_path}')

		# also write a short plain CSV with top-N anomalies
		topn = max(1, int(math.ceil(0.05 * combined.shape[0])))
		topn_path = os.path.join(args.out_dir, 'top_anomalies.csv')
		combined.head(topn).to_csv(topn_path)
		print(f'Wrote top {topn} anomalies to {topn_path}')

		# produce richer explanations and enrichment analysis
		explain_outdir = args.explain_output if args.explain_output else args.out_dir
		expl_path = None
		try:
			expl_path = write_explanations(df, combined, explain_outdir, top_n=args.explain_top_n)
			print(f'Wrote per-sample explanations to {expl_path}')
		except Exception as e:
			print(f'Failed to write explanations: {e}', file=sys.stderr)
		# compute enrichment across all samples and write CSV
		try:
			enr = compute_enrichment({s: {int(h): int(v) for h, v in df[s].fillna(0).astype(int).to_dict().items()} for s in df.columns})
			enr_path = os.path.join(explain_outdir, 'hash_enrichment.csv')
			enr.to_csv(enr_path)
			print(f'Wrote enrichment results to {enr_path}')
		except Exception as e:
			print(f'Failed enrichment computation: {e}', file=sys.stderr)

	elif mode_used == 'report':
		# Use the OutliMer report to produce heuristic explanations and rankings
		explained = explain_report_based(report_df)
		out_report = os.path.join(args.out_dir, 'outliers_report_from_reportcsv.csv')
		explained.to_csv(out_report)
		print(f'Wrote outlier report (from OutliMer report) to {out_report}')
		# also write a human-readable summary
		text_path = os.path.join(args.out_dir, 'outlier_explanations.txt')
		with open(text_path, 'w') as tf:
			for s in explained.index:
				row = explained.loc[s]
				tf.write(f"{s}: score={row['reason_score']:.4f}; reasons={row['reasons']}; pct_in_db={row['pct_in_db']:.3f}; n_hashes={int(row['n_hashes'])}; n_new_hashes={int(row['n_new_hashes'])}\n")
		print(f'Wrote textual explanations to {text_path}')

	else:
		raise RuntimeError('No valid mode selected')

	return 0


if __name__ == '__main__':
	raise SystemExit(main())



