"""Core OutliMer sketching and reporting command-line interface."""

import argparse
import concurrent.futures
import csv
import gzip
import importlib.metadata
import json
import logging
import os
import platform
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from typing import Dict, Iterator, List, Optional, Set, Tuple

from OutliMer import __version__

SOURMASH_IMPORT_ERROR: Optional[Exception] = None
try:
    import sourmash
    from sourmash import MinHash
except Exception as exc:
    sourmash = None
    MinHash = None
    SOURMASH_IMPORT_ERROR = exc

# ─────────────────────────────────────────────────────────────────────────────
# Module-level constants  (R2: replaces duplicated inline dicts)
# ─────────────────────────────────────────────────────────────────────────────

FASTQ_EXTENSIONS: List[str] = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
FASTA_EXTENSIONS: List[str] = [
    '.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz',
]
SIG_EXTENSIONS: List[str] = ['.sig', '.sig.gz']

EXTENSIONS_MAP: Dict[str, List[str]] = {
    'paired-fastq': FASTQ_EXTENSIONS,
    'single-fastq': FASTQ_EXTENSIONS,
    'fasta': FASTA_EXTENSIONS,
    'signature': SIG_EXTENSIONS,
}

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def open_maybe_gz(path: str):
    """Open a plain or gzip-compressed text file for reading."""
    if path.lower().endswith('.gz'):
        return gzip.open(path, 'rt', errors='strict')
    return open(path, 'r', errors='strict')


def fastq_sequences(path: str) -> Iterator[str]:
    """Yield sequences from a validated four-line FASTQ file."""
    with open_maybe_gz(path) as fh:
        record = 0
        while True:
            header = fh.readline()
            if not header:
                break
            record += 1
            seq = fh.readline()
            plus = fh.readline()
            quality = fh.readline()
            if not seq or not plus or not quality:
                raise ValueError(
                    f'{path}: truncated FASTQ record {record}')
            if not header.startswith('@'):
                raise ValueError(
                    f'{path}: record {record} header does not start with @')
            if not plus.startswith('+'):
                raise ValueError(
                    f'{path}: record {record} separator does not start with +')
            sequence = seq.rstrip('\r\n')
            qualities = quality.rstrip('\r\n')
            if not sequence:
                raise ValueError(f'{path}: record {record} has an empty sequence')
            if len(sequence) != len(qualities):
                raise ValueError(
                    f'{path}: record {record} sequence/quality lengths differ '
                    f'({len(sequence)} != {len(qualities)})')
            yield sequence


def fasta_sequences(path: str) -> Iterator[str]:
    """Yield sequences from a validated FASTA file."""
    with open_maybe_gz(path) as fh:
        seq_lines: List[str] = []
        seen_header = False
        record = 0
        for line in fh:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                if seen_header and not seq_lines:
                    raise ValueError(
                        f'{path}: FASTA record {record} has no sequence')
                if seen_header:
                    yield ''.join(seq_lines)
                record += 1
                seen_header = True
                seq_lines = []
                continue
            if not seen_header:
                raise ValueError(f'{path}: sequence found before first FASTA header')
            seq_lines.append(line.strip())
        if not seen_header:
            raise ValueError(f'{path}: no FASTA records found')
        if not seq_lines:
            raise ValueError(f'{path}: FASTA record {record} has no sequence')
        yield ''.join(seq_lines)


def load_sample_metadata(
    path: str,
    sample_column: str = 'sample',
) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """Load sample metadata from CSV/TSV keyed by sample name."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    delimiter = '\t' if path.lower().endswith(('.tsv', '.tab')) else ','
    with open(path, newline='') as fh:
        sample_text = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample_text, delimiters=',\t')
            delimiter = dialect.delimiter
        except csv.Error:
            pass
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError('metadata file has no header row')
        if sample_column not in reader.fieldnames:
            raise ValueError(
                f'metadata sample column {sample_column!r} not found')
        metadata_columns = [c for c in reader.fieldnames if c != sample_column]
        metadata: Dict[str, Dict[str, str]] = {}
        for row in reader:
            sample = (row.get(sample_column) or '').strip()
            if not sample:
                raise ValueError('metadata contains a blank sample name')
            if sample in metadata:
                raise ValueError(f'duplicate metadata sample: {sample}')
            metadata[sample] = {
                c: (row.get(c) or '') for c in metadata_columns
            }
    return metadata, metadata_columns


# ─────────────────────────────────────────────────────────────────────────────
# K-mer utilities
# ─────────────────────────────────────────────────────────────────────────────

def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(trans)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of kmer and its reverse complement."""
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


# ─────────────────────────────────────────────────────────────────────────────
# MinHash helpers
# ─────────────────────────────────────────────────────────────────────────────

def _extract_hash_counts(mh: 'MinHash') -> Dict[int, int]:
    """Extract {hash: abundance} from a sourmash MinHash object.

    R1: consolidated from duplicated logic in build_minhash_for_pair and
    build_minhash_for_single.
    B2: the dead second `if hasattr(mh, "hashes")` fallback is removed;
    if the first attempt raises, we try get_mins/get_abundance, then raise.
    """
    if hasattr(mh, 'hashes'):
        try:
            return {int(k): int(v) for k, v in mh.hashes.items()}
        except Exception as exc:
            logging.getLogger('OutliMer').debug(
                'mh.hashes dict access failed: %s', exc)

    if hasattr(mh, 'get_mins') and hasattr(mh, 'get_abundance'):
        try:
            return {int(h): int(mh.get_abundance(h)) for h in mh.get_mins()}
        except Exception as exc:
            logging.getLogger('OutliMer').debug(
                'get_mins/get_abundance failed: %s', exc)

    raise RuntimeError(
        'Could not extract hash counts from MinHash object. '
        'Check sourmash version compatibility.')


def _load_signature_counts(
    path: str,
    ksize: int,
    scaled: int,
    seed: int,
) -> Dict[int, int]:
    """Load one compatible DNA abundance signature at a common scale."""
    if sourmash is None:
        raise RuntimeError('sourmash not installed')

    if hasattr(sourmash, 'load_file_as_signatures'):
        try:
            signatures = list(sourmash.load_file_as_signatures(
                path, ksize=ksize, select_moltype='DNA'))
        except TypeError:
            signatures = list(sourmash.load_file_as_signatures(
                path, ksize=ksize))
    else:
        from sourmash.signature import load_signatures_from_json
        with open_maybe_gz(path) as fh:
            signatures = list(load_signatures_from_json(fh.read(), ksize=ksize))

    compatible = []
    for sig in signatures:
        mh = getattr(sig, 'minhash', None)
        if mh is None:
            continue
        sig_ksize = getattr(mh, 'ksize', ksize)
        if int(sig_ksize) != int(ksize):
            continue
        moltype = str(getattr(mh, 'moltype', 'DNA')).upper()
        if moltype != 'DNA':
            continue
        compatible.append(mh)

    if not compatible:
        raise RuntimeError(f'No DNA signature with ksize={ksize} found in {path}')
    if len(compatible) > 1:
        raise RuntimeError(
            f'{path} contains {len(compatible)} compatible signatures; '
            'provide one DNA signature per file')

    mh = compatible[0]
    sig_seed = int(getattr(mh, 'seed', seed))
    if sig_seed != seed:
        raise RuntimeError(
            f'{path} uses seed={sig_seed}, expected seed={seed}')
    if not bool(getattr(mh, 'track_abundance', True)):
        raise RuntimeError(
            f'{path} does not track hash abundance; abundance signatures are required')
    sig_scaled = int(getattr(mh, 'scaled', scaled))
    if sig_scaled <= 0:
        raise RuntimeError(f'{path} is not a scaled MinHash signature')
    if sig_scaled > scaled:
        raise RuntimeError(
            f'{path} uses scaled={sig_scaled}, which is sparser than requested '
            f'scaled={scaled} and cannot be upsampled')
    if sig_scaled < scaled:
        if not hasattr(mh, 'downsample'):
            raise RuntimeError(
                f'{path} uses scaled={sig_scaled}; sourmash cannot downsample '
                f'it to requested scaled={scaled}')
        mh = mh.downsample(scaled=scaled)

    return _extract_hash_counts(mh)


def _sketch_sample(
    name: str,
    r1: str,
    r2: str,
    ksize: int,
    scaled: int,
    seed: int,
    input_type: str,
) -> Tuple[str, Optional[Dict[int, int]], Optional[str]]:
    """Build a MinHash sketch for one sample.  Returns (name, counts, error).

    Designed to run inside a ThreadPoolExecutor (O4).
    sourmash releases the GIL during its C-level murmur hashing so threads
    provide genuine concurrency here.
    """
    if sourmash is None:
        return name, None, 'sourmash not installed'
    try:
        if input_type == 'signature':
            return name, _load_signature_counts(r1, ksize, scaled, seed), None

        mh = MinHash(
            ksize=ksize,
            n=0,
            scaled=scaled,
            seed=seed,
            track_abundance=True,
        )
        if r2:
            for seq in fastq_sequences(r1):
                mh.add_sequence(seq, force=True)
            for seq in fastq_sequences(r2):
                mh.add_sequence(seq, force=True)
        else:
            seq_iter = fasta_sequences if input_type == 'fasta' else fastq_sequences
            for seq in seq_iter(r1):
                mh.add_sequence(seq, force=True)
        return name, _extract_hash_counts(mh), None
    except Exception as exc:
        return name, None, str(exc)


# ─────────────────────────────────────────────────────────────────────────────
# Sample discovery
# ─────────────────────────────────────────────────────────────────────────────

def gather_samples_from_dir(
    input_dir: str,
    input_type: str = 'paired-fastq',
    substring: Optional[str] = None,
) -> List[Tuple[str, str, str]]:
    """Scan *input_dir* and return [(sample_name, r1, r2)] tuples.

    R2: uses the module-level EXTENSIONS_MAP instead of a locally-redefined dict.
    For single-fastq and fasta types, r2 is always ''.
    """
    if not os.path.isdir(input_dir):
        raise ValueError(f'input directory does not exist: {input_dir}')
    allowed_exts = EXTENSIONS_MAP.get(input_type, FASTQ_EXTENSIONS)

    entries: List[str] = []
    for fn in sorted(os.listdir(input_dir)):
        if substring and substring not in fn:
            continue
        fp = os.path.join(input_dir, fn)
        if not os.path.isfile(fp):
            continue
        if any(fn.lower().endswith(ext) for ext in allowed_exts):
            entries.append(fn)

    if input_type == 'paired-fastq':
        pairs_map: Dict[str, Dict[str, object]] = {}
        for fn in entries:
            stem = fn
            for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
                if stem.lower().endswith(ext):
                    stem = stem[:-len(ext)]
                    break
            read: Optional[str] = None
            prefix: Optional[str] = None
            tokens = ['_R1', '_R2', '_r1', '_r2', '_1', '_2',
                      '.R1', '.R2', '.r1', '.r2']
            for tok in tokens:
                pos = stem.rfind(tok)
                if pos != -1:
                    prefix = stem[:pos]
                    read = '1' if '1' in tok else '2'
                    break
            if prefix is None:
                m = re.match(r'(?i)(?P<prefix>.*?)[._-](?P<read>[12])(?:$|[._-])',
                             stem)
                if m:
                    prefix = m.group('prefix')
                    read = m.group('read')
            if prefix is None:
                prefix = stem
            if not prefix:
                raise ValueError(f'could not derive a sample name from {fn}')
            fullpath = os.path.join(input_dir, fn)
            rec = pairs_map.setdefault(prefix, {})
            if read == '1':
                if '1' in rec:
                    raise ValueError(
                        f'duplicate R1 files for sample {prefix}: '
                        f'{os.path.basename(str(rec["1"]))}, {fn}')
                rec['1'] = fullpath
            elif read == '2':
                if '2' in rec:
                    raise ValueError(
                        f'duplicate R2 files for sample {prefix}: '
                        f'{os.path.basename(str(rec["2"]))}, {fn}')
                rec['2'] = fullpath
            else:
                singles = rec.setdefault('singles', [])
                assert isinstance(singles, list)
                singles.append(fullpath)

        out: List[Tuple[str, str, str]] = []
        problems: List[str] = []
        for pfx, rec in sorted(pairs_map.items()):
            if '1' in rec and '2' in rec:
                out.append((pfx, str(rec['1']), str(rec['2'])))
            else:
                if '1' in rec or '2' in rec:
                    present = 'R1' if '1' in rec else 'R2'
                    problems.append(f'{pfx} has {present} but no matching mate')
                for path in rec.get('singles', []):
                    problems.append(
                        f'cannot determine R1/R2 from {os.path.basename(str(path))}')
        if problems:
            raise ValueError('invalid paired FASTQ inputs: ' + '; '.join(problems))
        return out

    # single-fastq or fasta
    out = []
    for fn in entries:
        stem = fn
        for ext in allowed_exts:
            if stem.lower().endswith(ext):
                stem = stem[:-len(ext)]
                break
        out.append((stem, os.path.join(input_dir, fn), ''))
    return out


def validate_sample_inputs(
    samples: List[Tuple[str, str, str]],
    input_type: str,
) -> None:
    """Validate sample names, file paths, and paired-read completeness."""
    seen: Set[str] = set()
    for name, r1, r2 in samples:
        if not name.strip():
            raise ValueError('sample names must not be blank')
        if name in seen:
            raise ValueError(f'duplicate sample name: {name}')
        seen.add(name)
        if not os.path.isfile(r1):
            raise ValueError(f'{name}: input file does not exist: {r1}')
        if input_type == 'paired-fastq':
            if not r2:
                raise ValueError(f'{name}: paired FASTQ sample has no R2 file')
            if not os.path.isfile(r2):
                raise ValueError(f'{name}: input file does not exist: {r2}')
            if os.path.abspath(r1) == os.path.abspath(r2):
                raise ValueError(f'{name}: R1 and R2 refer to the same file')
        elif r2:
            raise ValueError(
                f'{name}: R2 is only valid with --input-type paired-fastq')


# ─────────────────────────────────────────────────────────────────────────────
# CSV writers
# ─────────────────────────────────────────────────────────────────────────────

def write_wide_csv(
    sample_hash_counts: Dict[str, Dict[int, int]],
    out_path: str,
) -> None:
    all_hashes = sorted({h for c in sample_hash_counts.values() for h in c})
    samples = sorted(sample_hash_counts)
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['hash'] + samples)
        for h in all_hashes:
            writer.writerow([h] + [sample_hash_counts[s].get(h, 0) for s in samples])


def write_wide_csv_with_seq(
    sample_hash_counts: Dict[str, Dict[int, int]],
    out_path: str,
    seq_map: Dict[int, List[Tuple[str, str]]],
) -> None:
    all_hashes = sorted({h for c in sample_hash_counts.values() for h in c})
    samples = sorted(sample_hash_counts)
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['hash', 'sequence'] + samples)
        for h in all_hashes:
            seqs = seq_map.get(h, [])
            seq_field = ';'.join(k for _, k in seqs)
            writer.writerow([h, seq_field] + [sample_hash_counts[s].get(h, 0)
                                               for s in samples])


def write_long_csv(
    sample_hash_counts: Dict[str, Dict[int, int]],
    out_path: str,
) -> None:
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'hash', 'count'])
        for sample in sorted(sample_hash_counts):
            for h, c in sample_hash_counts[sample].items():
                if c:
                    writer.writerow([sample, h, c])


def write_long_csv_with_seq(
    sample_hash_counts: Dict[str, Dict[int, int]],
    out_path: str,
    seq_map: Dict[int, List[Tuple[str, str]]],
) -> None:
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'hash', 'sequence', 'count'])
        for sample in sorted(sample_hash_counts):
            for h, c in sample_hash_counts[sample].items():
                if c:
                    seqs = seq_map.get(h, [])
                    seq_field = ';'.join(k for _, k in seqs)
                    writer.writerow([sample, h, seq_field, c])


def write_union_summary(
    sample_hash_counts: Dict[str, Dict[int, int]],
    out_dir: Optional[str],
    seq_map: Optional[Dict[int, List[Tuple[str, str]]]] = None,
    filename: str = 'top_union_summary.csv',
    out_path: Optional[str] = None,
) -> str:
    """Write a union summary CSV of all hashes across all samples."""
    if not sample_hash_counts:
        raise RuntimeError('No sample hash counts available')
    if out_path is None:
        export_dir = os.path.join(out_dir or '.', 'export_kmers')
        os.makedirs(export_dir, exist_ok=True)
        out_path = os.path.join(export_dir, filename)
    else:
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    all_hashes = sorted({int(h) for c in sample_hash_counts.values()
                          for h in c})
    samples = sorted(sample_hash_counts)
    include_seq = bool(seq_map)
    header = ['hash'] + (['sequence'] if include_seq else []) + samples
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for h in all_hashes:
            row: list = [h]
            if include_seq:
                seqs = seq_map.get(h, [])  # type: ignore[union-attr]
                row.append(seqs[0][1] if seqs else '')
            for s in samples:
                row.append(int(sample_hash_counts[s].get(h, 0)))
            writer.writerow(row)
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# Hash database
# ─────────────────────────────────────────────────────────────────────────────

def load_hash_db(path: str) -> Set[int]:
    result: Set[int] = set()
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                result.add(int(line))
            except ValueError:
                continue
    return result


def save_hash_db(path: str, hashes: Set[int]) -> None:
    open_fn = gzip.open if path.endswith('.gz') else open
    mode = 'wt' if path.endswith('.gz') else 'w'
    with open_fn(path, mode) as fh:  # type: ignore[call-overload]
        for h in sorted(hashes):
            fh.write(str(h) + '\n')


# ─────────────────────────────────────────────────────────────────────────────
# Comparison
# ─────────────────────────────────────────────────────────────────────────────

def compare_sample_to_db(
    sample_hashes: Set[int],
    db_hashes: Set[int],
) -> Tuple[int, int, float]:
    """Return (n_in_db, n_total, fraction_in_db).

    O1: uses set intersection instead of an explicit loop.
    """
    n_total = len(sample_hashes)
    if n_total == 0:
        return 0, 0, 0.0
    n_in = len(sample_hashes & db_hashes)
    return n_in, n_total, n_in / n_total


# ─────────────────────────────────────────────────────────────────────────────
# K-mer example finder
# ─────────────────────────────────────────────────────────────────────────────

def find_examples_for_hashes(
    sample: str,
    r1: str,
    r2: str,
    target_hashes: Set[int],
    ksize: int,
    per_hash: int = 1,
) -> Dict[int, List[Tuple[str, str]]]:
    """Scan paired FASTQ for canonical k-mers that map to *target_hashes*.

    B5: guarded against empty r2.
    O3: inner k-mer loop breaks early once all target hashes are satisfied.
    R5: seq.strip() removed – fastq_sequences already strips.
    """
    out: Dict[int, List[Tuple[str, str]]] = {}
    if not target_hashes or not r1:
        return out
    try:
        import sourmash.minhash as sm
        seed = sm.get_minhash_default_seed()
        max_hash = sm.get_minhash_max_hash()
    except Exception as exc:
        logging.getLogger('OutliMer').debug(
            'Could not import sourmash.minhash: %s', exc)
        return out

    def _h(kmer: str) -> int:
        return int(sm.hash_murmur(kmer, seed)) & max_hash

    remaining = set(target_hashes)

    def _scan(seq_iter: Iterator[str]) -> None:
        for seq in seq_iter:
            L = len(seq)
            if L < ksize:
                continue
            for i in range(L - ksize + 1):
                kmer = seq[i:i + ksize]
                if 'N' in kmer or 'n' in kmer:
                    continue
                kcanon = canonical_kmer(kmer)
                h = _h(kcanon)
                if h in target_hashes:
                    lst = out.setdefault(h, [])
                    if len(lst) < per_hash:
                        lst.append((sample, kcanon))
                        if len(lst) >= per_hash:
                            remaining.discard(h)
                # O3: break inner loop as soon as all hashes are found
                if not remaining:
                    return
            if not remaining:
                return

    _scan(fastq_sequences(r1))
    if remaining and r2:  # B5: only scan r2 when it is a real path
        _scan(fastq_sequences(r2))

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Shared computation helpers
# ─────────────────────────────────────────────────────────────────────────────

def _build_kmer_to_samples(
    sample_hash_counts: Dict[str, Dict[int, int]],
    exported_items_per_sample: Dict[str, List[Tuple[int, int]]],
    export_mode: str,
) -> Dict[int, Dict[str, int]]:
    """Build hash -> {sample: count} mapping.

    R6: centralised so the same aggregation is not repeated in multiple places.
    """
    kmer_to_samples: Dict[int, Dict[str, int]] = defaultdict(dict)
    if export_mode == 'full':
        for sample, counts in sample_hash_counts.items():
            for h, c in counts.items():
                kmer_to_samples[int(h)][sample] = int(c)
    else:
        for sample, items in exported_items_per_sample.items():
            for h, c in items:
                kmer_to_samples[int(h)][sample] = int(c)
    return kmer_to_samples


def _compute_total_counts(
    sample_hash_counts: Dict[str, Dict[int, int]],
) -> Dict[int, int]:
    """Sum hash abundances across all samples.

    R6: computed once and reused wherever a global total is needed.
    """
    total: Dict[int, int] = {}
    for counts in sample_hash_counts.values():
        for h, c in counts.items():
            total[int(h)] = total.get(int(h), 0) + int(c)
    return total


def _safe_cache_name(sample: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', sample).strip('._') or 'sample'


def _source_fingerprint(paths: Tuple[str, str]) -> List[Dict[str, int | str]]:
    fingerprint: List[Dict[str, int | str]] = []
    for path in paths:
        if not path:
            continue
        stat = os.stat(path)
        fingerprint.append({
            'path': os.path.abspath(path),
            'size': int(stat.st_size),
            'mtime_ns': int(stat.st_mtime_ns),
        })
    return fingerprint


def _runtime_versions() -> Dict[str, str]:
    versions = {'OutliMer': __version__, 'Python': platform.python_version()}
    for distribution in (
        'sourmash', 'numpy', 'pandas', 'scipy', 'scikit-learn',
        'matplotlib', 'seaborn',
    ):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = 'not installed'
    return versions


def _cache_payload_metadata(
    sample: str,
    paths: Tuple[str, str],
    ksize: int,
    scaled: int,
    seed: int,
    input_type: str,
) -> Dict[str, object]:
    return {
        'sample': sample,
        'inputs': _source_fingerprint(paths),
        'ksize': int(ksize),
        'scaled': int(scaled),
        'seed': int(seed),
        'input_type': input_type,
        'schema': 2,
    }


def _cache_path(cache_dir: str, sample: str) -> str:
    return os.path.join(cache_dir, f'{_safe_cache_name(sample)}.json')


def _load_cached_counts(
    cache_dir: Optional[str],
    sample: str,
    paths: Tuple[str, str],
    ksize: int,
    scaled: int,
    seed: int,
    input_type: str,
) -> Optional[Dict[int, int]]:
    if not cache_dir:
        return None
    path = _cache_path(cache_dir, sample)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        payload = json.load(fh)
    expected = _cache_payload_metadata(
        sample, paths, ksize, scaled, seed, input_type)
    if payload.get('metadata') != expected:
        return None
    counts = payload.get('counts', {})
    return {int(h): int(c) for h, c in counts.items()}


def _save_cached_counts(
    cache_dir: Optional[str],
    sample: str,
    paths: Tuple[str, str],
    ksize: int,
    scaled: int,
    seed: int,
    input_type: str,
    counts: Dict[int, int],
) -> None:
    if not cache_dir:
        return
    os.makedirs(cache_dir, exist_ok=True)
    payload = {
        'metadata': _cache_payload_metadata(
            sample, paths, ksize, scaled, seed, input_type),
        'counts': {str(int(h)): int(c) for h, c in counts.items()},
    }
    with open(_cache_path(cache_dir, sample), 'w') as fh:
        json.dump(payload, fh, sort_keys=True)


def write_multiqc_summary(
    report_rows: List[Tuple[str, int, int, float, int]],
    out_dir: str,
) -> str:
    """Write MultiQC custom-content JSON for core OutliMer run stats."""
    payload = {
        'id': 'outlimer_sample_summary',
        'section_name': 'OutliMer sample summary',
        'description': 'OutliMer per-sample hash sharing and novelty metrics.',
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'outlimer_new_hashes',
            'title': 'OutliMer novel hashes',
            'ylab': 'Novel hashes',
        },
        'data': {
            sample: {
                'n_hashes': int(n_hashes),
                'n_in_db': int(n_in),
                'pct_in_db': float(pct),
                'n_new_hashes': int(n_new),
            }
            for sample, n_hashes, n_in, pct, n_new in report_rows
        },
    }
    out_path = os.path.join(out_dir, 'outlimer_summary_mqc.json')
    with open(out_path, 'w') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return out_path


def write_ro_crate_metadata(
    out_dir: str,
    inputs: Dict[str, str],
    outputs: Dict[str, Optional[str]],
    parameters: Dict[str, object],
) -> str:
    """Write minimal RO-Crate metadata for an OutliMer run."""
    file_entities = []
    for label, path in {**inputs, **outputs}.items():
        if not path:
            continue
        file_entities.append({
            '@id': os.path.basename(path),
            '@type': 'File',
            'name': label,
            'contentUrl': path,
        })
    crate = {
        '@context': 'https://w3id.org/ro/crate/1.1/context',
        '@graph': [
            {
                '@id': 'ro-crate-metadata.json',
                '@type': 'CreativeWork',
                'conformsTo': {'@id': 'https://w3id.org/ro/crate/1.1'},
                'about': {'@id': './'},
            },
            {
                '@id': './',
                '@type': 'Dataset',
                'name': 'OutliMer run',
                'dateCreated': datetime.now(timezone.utc).isoformat(),
                'hasPart': [{'@id': item['@id']} for item in file_entities],
                'softwareVersion': __version__,
                'outlimerParameters': parameters,
            },
            *file_entities,
        ],
    }
    out_path = os.path.join(out_dir, 'ro-crate-metadata.json')
    with open(out_path, 'w') as fh:
        json.dump(crate, fh, indent=2)
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# Sample processing (D2: extracted from main)
# ─────────────────────────────────────────────────────────────────────────────

def process_samples(
    pairs: List[Tuple[str, str, str]],
    ksize: int,
    scaled: int,
    seed: int,
    input_type: str,
    collect_sequences: bool,
    top_per_sample: int,
    profile_samples: Set[str],
    profile_dir: Optional[str],
    profile_force_recompute: bool,
    profile_agg_counts: Dict[int, int],
    profile_update: bool,
    profile_name: Optional[str],
    threads: int,
    cache_dir: Optional[str],
    logger: logging.Logger,
) -> Tuple[
    Dict[str, Dict[int, int]],   # sample_hash_counts
    Dict[str, Tuple[str, str]],  # sample_files
    Dict[int, List[Tuple[str, str]]],  # sourmash_seq_map
    Dict[int, int],              # updated profile_agg_counts
    Set[str],                    # updated profile_samples
    Dict[str, str],              # failed sample -> error
]:
    """Sketch all samples in parallel and collect sequence mappings.

    O4: uses ThreadPoolExecutor so multiple samples are sketched concurrently.
    D4: exception details are now logged rather than silently swallowed.
    """
    sample_hash_counts: Dict[str, Dict[int, int]] = {}
    sample_files: Dict[str, Tuple[str, str]] = {}
    sourmash_seq_map: Dict[int, List[Tuple[str, str]]] = {}
    failures: Dict[str, str] = {}

    # Determine which samples need sketching
    to_sketch: List[Tuple[str, str, str]] = []
    for name, r1, r2 in pairs:
        if (profile_dir and name in profile_samples
                and not profile_force_recompute):
            logger.info('Skipping %s (already in profile)', name)
            continue
        try:
            cached = _load_cached_counts(
                cache_dir, name, (r1, r2), ksize, scaled, seed, input_type)
        except Exception as exc:
            logger.debug('Cache read failed for %s: %s', name, exc)
            cached = None
        if cached is not None:
            sample_hash_counts[name] = cached
            sample_files[name] = (r1, r2)
            total_count = sum(cached.values())
            logger.info('CACHE %s: unique_hashes=%d total_count=%d',
                        name, len(cached), total_count)
            if profile_update and profile_name:
                for h, c in cached.items():
                    profile_agg_counts[int(h)] = (
                        profile_agg_counts.get(int(h), 0) + int(c))
                profile_samples.add(name)
            if collect_sequences and cached and input_type != 'signature':
                top_n = top_per_sample if top_per_sample > 0 else None
                items = sorted(cached.items(), key=lambda x: -x[1])
                targets = {h for h, _ in (items[:top_n] if top_n else items)}
                if targets:
                    try:
                        found = find_examples_for_hashes(
                            name, r1, r2, targets, ksize, per_hash=1)
                        for h, lst in found.items():
                            cur = sourmash_seq_map.setdefault(int(h), [])
                            if not cur and lst:
                                cur.append(lst[0])
                    except Exception as exc:
                        logger.debug('Sequence mapping failed for %s: %s',
                                     name, exc)
            continue
        to_sketch.append((name, r1, r2))

    logger.info('Sketching %d sample(s) with %d thread(s) …', len(to_sketch),
                threads)

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(
                _sketch_sample, name, r1, r2, ksize, scaled, seed, input_type):
                (name, r1, r2)
            for name, r1, r2 in to_sketch
        }
        for fut in concurrent.futures.as_completed(futures):
            name, r1, r2 = futures[fut]
            try:
                sketch_name, counts, error = fut.result()
            except Exception as exc:
                sketch_name, counts, error = name, None, str(exc)
            if error or counts is None:
                logger.warning('Failed to process sample %s: %s', name, error)
                failures[name] = error or 'unknown processing error'
                continue
            sample_hash_counts[name] = counts
            sample_files[name] = (r1, r2)
            try:
                _save_cached_counts(
                    cache_dir, name, (r1, r2), ksize, scaled, seed,
                    input_type, counts)
            except Exception as exc:
                logger.debug('Cache write failed for %s: %s', name, exc)

            total_count = sum(counts.values())
            logger.info('DONE %s: unique_hashes=%d total_count=%d',
                        name, len(counts), total_count)

            # Profile matching stats
            if profile_agg_counts:
                n_in = sum(1 for h in counts if int(h) in profile_agg_counts)
                frac = n_in / len(counts) if counts else 0.0
                logger.info('  %s: %d/%d hashes in profile (frac=%.4f)',
                            name, n_in, len(counts), frac)

            # Profile update
            if profile_update and profile_name:
                for h, c in counts.items():
                    profile_agg_counts[int(h)] = (
                        profile_agg_counts.get(int(h), 0) + int(c))
                profile_samples.add(name)

            # Optional per-sample sequence collection
            if collect_sequences and counts:
                top_n = top_per_sample if top_per_sample > 0 else None
                items = sorted(counts.items(), key=lambda x: -x[1])
                targets = {h for h, _ in (items[:top_n] if top_n else items)}
                if targets:
                    try:
                        found = find_examples_for_hashes(
                            name, r1, r2, targets, ksize, per_hash=1)
                        for h, lst in found.items():
                            cur = sourmash_seq_map.setdefault(int(h), [])
                            if not cur and lst:
                                cur.append(lst[0])
                    except Exception as exc:
                        logger.debug('Sequence mapping failed for %s: %s',
                                     name, exc)

    return (sample_hash_counts, sample_files, sourmash_seq_map,
            profile_agg_counts, profile_samples, failures)


# ─────────────────────────────────────────────────────────────────────────────
# Per-sample file exports (D2: extracted from main)
# ─────────────────────────────────────────────────────────────────────────────

def _write_per_sample_top_n(
    sample: str,
    items: List[Tuple[int, int]],
    kmer_to_samples: Dict[int, Dict[str, int]],
    sourmash_seq_map: Dict[int, List[Tuple[str, str]]],
    seq_enabled: bool,
    outdir: str,
    top_n: Optional[int],
    fill_unique: bool,
    write_unique_shared: bool,
    logger: logging.Logger,
) -> None:
    """Write {sample}_topN.csv, _unique_hashes.csv, _shared_hashes.csv."""
    unique_rows: List[tuple] = []
    shared_rows: List[tuple] = []
    for h, cnt in items:
        owners = kmer_to_samples.get(int(h), {})
        seq = ''
        if seq_enabled:
            seqs = sourmash_seq_map.get(int(h), [])
            seq = seqs[0][1] if seqs else ''
        is_unique = len(owners) == 1 and sample in owners
        other = ';'.join(f'{s}:{c}' for s, c in sorted(owners.items())
                         if s != sample)
        if is_unique:
            unique_rows.append((int(h), seq, int(cnt)))
        else:
            shared_rows.append((int(h), seq, int(cnt), other))

    # topN file
    samp_topN = os.path.join(outdir, f'{sample}_topN.csv')
    with open(samp_topN, 'w', newline='') as fh:
        writer = csv.writer(fh)
        if seq_enabled:
            writer.writerow(['hash', 'sequence', 'count', 'is_unique',
                             'other_samples_counts'])
            for h, seq, cnt in unique_rows:
                owners = kmer_to_samples.get(h, {})
                other = ';'.join(f'{s}:{c}' for s, c in sorted(owners.items())
                                 if s != sample)
                writer.writerow([h, seq, cnt, 'yes', ''])
            for h, seq, cnt, other in shared_rows:
                writer.writerow([h, seq, cnt, 'no', other])
        else:
            writer.writerow(['hash', 'count', 'is_unique',
                             'other_samples_counts'])
            for h, seq, cnt in unique_rows:
                owners = kmer_to_samples.get(h, {})
                other = ';'.join(f'{s}:{c}' for s, c in sorted(owners.items())
                                 if s != sample)
                writer.writerow([h, cnt, 'yes', ''])
            for h, seq, cnt, other in shared_rows:
                writer.writerow([h, cnt, 'no', other])
    logger.debug('Wrote %s (%d rows)', samp_topN, len(items))

    if not write_unique_shared:
        return

    # unique / shared files
    unique_file = os.path.join(outdir, f'{sample}_unique_hashes.csv')
    shared_file = os.path.join(outdir, f'{sample}_shared_hashes.csv')
    with open(unique_file, 'w', newline='') as uf, \
         open(shared_file, 'w', newline='') as sf:
        uw = csv.writer(uf)
        sw = csv.writer(sf)
        if seq_enabled:
            uw.writerow(['hash', 'sequence', 'count'])
            sw.writerow(['hash', 'sequence', 'count', 'other_samples_counts'])
        else:
            uw.writerow(['hash', 'count'])
            sw.writerow(['hash', 'count', 'other_samples_counts'])
        n_u = n_s = 0
        for h, seq, cnt in unique_rows:
            if top_n is None or n_u < top_n:
                uw.writerow([h, seq, cnt] if seq_enabled else [h, cnt])
                n_u += 1
        for h, seq, cnt, other in shared_rows:
            if top_n is None or n_s < top_n:
                sw.writerow([h, seq, cnt, other] if seq_enabled
                            else [h, cnt, other])
                n_s += 1
    logger.debug('Wrote %s (unique=%d) %s (shared=%d)',
                 unique_file, n_u, shared_file, n_s)

    # optional filled-unique file
    if fill_unique and top_n is not None:
        filled_file = os.path.join(outdir, f'{sample}_unique_hashes_filled.csv')
        with open(filled_file, 'w', newline='') as ff:
            fw = csv.writer(ff)
            if seq_enabled:
                fw.writerow(['hash', 'sequence', 'count', 'is_unique',
                             'other_samples_counts'])
            else:
                fw.writerow(['hash', 'count', 'is_unique',
                             'other_samples_counts'])
            written = 0
            for h, seq, cnt in unique_rows:
                if written >= top_n:
                    break
                fw.writerow([h, seq, cnt, 'yes', ''] if seq_enabled
                            else [h, cnt, 'yes', ''])
                written += 1
            for h, seq, cnt, other in shared_rows:
                if written >= top_n:
                    break
                fw.writerow([h, seq, cnt, 'no', other] if seq_enabled
                            else [h, cnt, 'no', other])
                written += 1
        logger.debug('Wrote filled unique file %s (%d rows)', filled_file,
                     written)


def _write_shared_summary(
    exported_items: Dict[str, List[Tuple[int, int]]],
    kmer_to_samples: Dict[int, Dict[str, int]],
    out_path: str,
    top_shared: int,
) -> None:
    """Write top shared hashes for each sample."""
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    limit = top_shared if top_shared > 0 else None
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([
            'sample', 'hash', 'count', 'other_samples_counts',
            'n_other_samples',
        ])
        for sample in sorted(exported_items):
            shared: List[Tuple[int, int, str, int]] = []
            for h, cnt in exported_items[sample]:
                owners = kmer_to_samples.get(int(h), {})
                others = {s: c for s, c in owners.items() if s != sample}
                if not others:
                    continue
                other_field = ';'.join(
                    f'{s}:{c}' for s, c in sorted(others.items()))
                shared.append((int(h), int(cnt), other_field, len(others)))
            shared.sort(key=lambda row: (-row[3], -row[1], row[0]))
            for h, cnt, other_field, n_other in (
                    shared[:limit] if limit is not None else shared):
                writer.writerow([sample, h, cnt, other_field, n_other])


def _run_exports(
    args: argparse.Namespace,
    sample_hash_counts: Dict[str, Dict[int, int]],
    sample_files: Dict[str, Tuple[str, str]],
    sourmash_seq_map: Dict[int, List[Tuple[str, str]]],
    top_n: Optional[int],
    logger: logging.Logger,
) -> None:
    """Write per-sample CSV exports.  D2: extracted from main."""
    outdir = args.export_kmers_dir
    os.makedirs(outdir, exist_ok=True)
    export_mode = args.export_mode
    seq_enabled = args.collect_sequences or args.export_sequences

    # Build top-N exported items per sample
    exported_items: Dict[str, List[Tuple[int, int]]] = {}
    for sample, counts in sample_hash_counts.items():
        items = sorted(counts.items(), key=lambda x: -x[1])
        if top_n and top_n > 0:
            if len(items) < top_n:
                if args.force_top:
                    raise RuntimeError(
                        f'Sample {sample} has only {len(items)} hashes; '
                        f'cannot satisfy --force-top {top_n}')
                logger.debug(
                    'Sample %s: MinHash has only %d hashes (requested %d)',
                    sample, len(items), top_n)
            items = items[:top_n]
        exported_items[sample] = items

    # B3 / R6: compute kmer_to_samples once here (was computed again per branch)
    kmer_to_samples = _build_kmer_to_samples(
        sample_hash_counts, exported_items, export_mode)

    # Optionally build sequence mappings in parallel
    if args.collect_sequences:
        logger.info('Building sequence examples for top-%s hashes …', top_n)
        samples_to_map = [s for s in exported_items if exported_items[s]]
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.map_threads) as pool:
            futures = {
                pool.submit(
                    find_examples_for_hashes, s,
                    sample_files[s][0], sample_files[s][1],
                    {h for h, _ in exported_items[s]},
                    args.ksize, 1,
                ): s
                for s in samples_to_map
            }
            for fut in concurrent.futures.as_completed(futures):
                s = futures[fut]
                try:
                    found = fut.result()
                except Exception as exc:
                    logger.warning('Mapping failed for %s: %s', s, exc)
                    continue
                for h, lst in found.items():
                    cur = sourmash_seq_map.setdefault(int(h), [])
                    if not cur and lst:
                        cur.append(lst[0])
                logger.debug('Mapped %s: %d examples', s, len(found))

    # Optional hash->sequence mapping file
    if args.map_out:
        with open(args.map_out, 'w', newline='') as mf:
            mw = csv.writer(mf)
            if seq_enabled:
                mw.writerow(['hash', 'sample', 'sequence'])
                for h in sorted(sourmash_seq_map):
                    for samp, seq in sourmash_seq_map[h]:
                        mw.writerow([h, samp, seq])
            else:
                logger.warning('--map-out specified but no sequences collected')
                mw.writerow(['hash', 'sample'])
                for h in sorted(sourmash_seq_map):
                    for samp, _ in sourmash_seq_map[h]:
                        mw.writerow([h, samp])
        logger.info('Wrote hash->sequence mapping to %s', args.map_out)

    fast_mode = not seq_enabled

    # ── Per-sample top-N / unique / shared files ──────────────────────────────
    for sample in sorted(exported_items):
        _write_per_sample_top_n(
            sample, exported_items[sample], kmer_to_samples,
            sourmash_seq_map, seq_enabled, outdir, top_n,
            args.fill_unique, args.export_unique_shared or args.fill_unique,
            logger)

    # ── Lightweight global shared summary (fast mode) ─────────────────────────
    global_shared_path = (args.export_global_shared
                          or os.path.join(outdir, 'shared_hashes_all.csv'))
    with open(global_shared_path, 'w', newline='') as gh:
        gw = csv.writer(gh)
        gw.writerow(['hash', 'samples_counts', 'n_samples'])
        n_shared = 0
        for h, smap in sorted(kmer_to_samples.items(),
                               key=lambda x: -len(x[1])):
            if len(smap) <= 1:
                continue
            n_shared += 1
            gw.writerow([h,
                         ';'.join(f'{s}:{c}' for s, c in sorted(smap.items())),
                         len(smap)])
    logger.info('Wrote global shared hashes (%d) to %s', n_shared,
                global_shared_path)

    if args.shared_summary:
        _write_shared_summary(
            exported_items, kmer_to_samples, args.shared_summary,
            args.top_shared)
        logger.info('Wrote shared summary to %s', args.shared_summary)

    if fast_mode:
        logger.info('Fast mode: skipping heavy sequence-export steps.')
        return

    # ── Heavy path: sequence exports ──────────────────────────────────────────
    if args.export_sequences:
        logger.info('Exporting canonical k-mer sequences …')
        for sample in sorted(sample_files):
            r1, r2 = sample_files[sample]
            if export_mode == 'full':
                target_hashes = set(int(h)
                                    for h in sample_hash_counts.get(sample, {}))
                counts_map = {int(h): int(c)
                              for h, c in sample_hash_counts.get(sample, {}).items()}
            else:
                exported = exported_items.get(sample, [])
                target_hashes = {int(h) for h, _ in exported}
                counts_map = {int(h): int(c) for h, c in exported}
            if not target_hashes:
                continue
            seq_for_hash: Dict[int, str] = {
                int(h): sourmash_seq_map[int(h)][0][1]
                for h in target_hashes
                if sourmash_seq_map.get(int(h))
            }
            remaining_hashes = {h for h in target_hashes
                                 if h not in seq_for_hash}
            if remaining_hashes:
                try:
                    found = find_examples_for_hashes(
                        sample, r1, r2, remaining_hashes, args.ksize,
                        per_hash=1)
                    for h, lst in found.items():
                        if lst:
                            seq_for_hash[int(h)] = lst[0][1]
                except Exception as exc:
                    logger.debug('Sequence scan failed for %s: %s', sample, exc)

            items_seq = sorted(
                ((h, counts_map.get(h, 0), seq_for_hash.get(h, ''))
                 for h in target_hashes),
                key=lambda x: -x[1])
            if top_n:
                items_seq = items_seq[:top_n]

            samp_csv = os.path.join(outdir, f'{sample}_kmers.csv')
            with open(samp_csv, 'w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow(['kmer', 'hash', 'count'])
                for h, cnt, seq in items_seq:
                    writer.writerow([seq, h, cnt])

            unique_k = os.path.join(outdir, f'{sample}_unique_kmers.csv')
            shared_k = os.path.join(outdir, f'{sample}_shared_kmers.csv')
            with open(unique_k, 'w', newline='') as uf, \
                 open(shared_k, 'w', newline='') as sf:
                uw = csv.writer(uf)
                sw = csv.writer(sf)
                uw.writerow(['kmer', 'hash', 'count'])
                sw.writerow(['kmer', 'hash', 'count', 'other_samples_counts'])
                nu = ns = 0
                for h, cnt, seq in items_seq:
                    owners = kmer_to_samples.get(h, {})
                    others = [s for s in sorted(owners) if s != sample]
                    if not others:
                        uw.writerow([seq, h, cnt])
                        nu += 1
                    else:
                        sw.writerow([seq, h, cnt,
                                     ';'.join(f'{s}:{owners[s]}' for s in others)])
                        ns += 1
            logger.info('Wrote sequence exports for %s (unique=%d shared=%d)',
                        sample, nu, ns)

    # B3: export_top_union written once, OUTSIDE any per-sample loop
    if args.export_top_union:
        union_path = (args.top_union_summary
                      or os.path.join(outdir, 'top_union_summary.csv'))
        # R6: reuse kmer_to_samples already built above
        total_counts = _compute_total_counts(sample_hash_counts)
        all_samples = sorted(exported_items)
        with open(union_path, 'w', newline='') as uf:
            uw = csv.writer(uf)
            header = ['hash']
            if seq_enabled:
                header.append('sequence')
            header += ['total_count'] + all_samples
            uw.writerow(header)
            for h in sorted(total_counts, key=lambda x: -total_counts[x]):
                row: list = [h]
                if seq_enabled:
                    seqs = sourmash_seq_map.get(h, [])
                    row.append(seqs[0][1] if seqs else '')
                row.append(total_counts[h])
                for s in all_samples:
                    row.append(sample_hash_counts[s].get(h, 0)
                               if s in sample_hash_counts else 0)
                uw.writerow(row)
        logger.info('Wrote union top-N summary to %s', union_path)


# ─────────────────────────────────────────────────────────────────────────────
# Outlier detection (D2: extracted from main)
# ─────────────────────────────────────────────────────────────────────────────

def _detect_outliers(
    args: argparse.Namespace,
    report_rows: List[Tuple[str, int, int, float, int]],
    sample_hash_counts: Dict[str, Dict[int, int]],
    sample_files: Dict[str, Tuple[str, str]],
    sourmash_seq_map: Dict[int, List[Tuple[str, str]]],
    db_hashes: Set[int],
    db_loaded: bool,
    out_dir: Optional[str],
    logger: logging.Logger,
) -> None:
    """Identify and export hashes unique to outlier samples.  D2."""
    threshold = args.outlier_threshold
    outliers = [r[0] for r in report_rows if float(r[3]) <= threshold]
    outlier_set = set(outliers)
    logger.info('Outlier samples (pct_in_db <= %.2f): %s', threshold, outliers)

    # hash -> set of samples containing it
    hash_to_samples: Dict[int, Set[str]] = defaultdict(set)
    for s, counts in sample_hash_counts.items():
        for h in counts:
            hash_to_samples[int(h)].add(s)

    if db_loaded:
        baseline = set(db_hashes)
    else:
        baseline = set()
        for s, counts in sample_hash_counts.items():
            if s not in outlier_set:
                baseline.update(counts)

    outlier_only = [h for h, samples in hash_to_samples.items()
                    if samples.issubset(outlier_set) and h not in baseline]

    # R6: top-N filter using pre-computed totals
    if args.outlier_top > 0:
        total_counts = _compute_total_counts(sample_hash_counts)
        outlier_only.sort(key=lambda x: -total_counts.get(x, 0))
        outlier_only = outlier_only[:args.outlier_top]

    outlier_dir = args.outlier_dir or os.path.join(out_dir or '.', 'outliers')
    os.makedirs(outlier_dir, exist_ok=True)

    hashes_file = os.path.join(outlier_dir, 'outlier_only_hashes.csv')
    with open(hashes_file, 'w', newline='') as hf:
        hw = csv.writer(hf)
        hw.writerow(['hash', 'samples_present'])
        for h in sorted(outlier_only):
            hw.writerow([h,
                         ';'.join(sorted(hash_to_samples.get(h, set())))])
    logger.info('Wrote outlier-only hashes to %s (%d)', hashes_file,
                len(outlier_only))

    for s in outliers:
        s_hashes = set(int(h) for h in sample_hash_counts.get(s, {}))
        new_vs_baseline = sorted(h for h in s_hashes if h not in baseline)
        samp_file = os.path.join(outlier_dir, f'{s}_outlier_new_hashes.csv')
        with open(samp_file, 'w', newline='') as sf:
            sw = csv.writer(sf)
            sw.writerow(['hash', 'count'])
            for h in new_vs_baseline:
                sw.writerow([h, int(sample_hash_counts[s].get(h, 0))])
        logger.info('Wrote outlier-new hashes for %s to %s (%d)',
                    s, samp_file, len(new_vs_baseline))

    if not args.outlier_seq_per_hash:
        return

    # Collect example sequences for outlier-only hashes
    seq_map: Dict[int, str] = {}
    for h in outlier_only:
        entries = sourmash_seq_map.get(h, [])
        if entries:
            seq_map[h] = entries[0][1]

    remaining_hashes = [h for h in outlier_only if h not in seq_map]
    if remaining_hashes:
        for s in outliers:
            r1, r2 = sample_files.get(s, ('', ''))
            if not r1:
                continue
            try:
                found = find_examples_for_hashes(
                    s, r1, r2, set(remaining_hashes), args.ksize, per_hash=1)
                for h, lst in found.items():
                    if lst:
                        seq_map[int(h)] = lst[0][1]
            except Exception as exc:
                logger.debug('Outlier seq scan failed for %s: %s', s, exc)

    seqs_file = os.path.join(outlier_dir, 'outlier_only_hashes_with_seq.csv')
    with open(seqs_file, 'w', newline='') as sf:
        sw = csv.writer(sf)
        sw.writerow(['hash', 'sequence', 'samples_present'])
        for h in sorted(outlier_only):
            sw.writerow([h, seq_map.get(h, ''),
                         ';'.join(sorted(hash_to_samples.get(h, set())))])
    logger.info('Wrote outlier hashes + sequences to %s', seqs_file)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting (D2: moved to module level from nested-inside-main)
# ─────────────────────────────────────────────────────────────────────────────

def _load_plotting():
    """Import plotting dependencies only when plots are requested."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt_mod
        import numpy as np_mod
        return plt_mod, np_mod
    except Exception as exc:
        logging.getLogger('OutliMer').warning(
            'matplotlib/numpy not available; cannot generate plots: %s', exc)
        return None, None


def generate_plots(
    report_rows: List[Tuple[str, int, int, float, int]],
    sample_hash_counts: Dict[str, Dict[int, int]],
    plot_prefix: str,
    db_loaded: bool,
    logger: logging.Logger,
) -> None:
    """Generate shared/unique stacked bar, Jaccard heatmap, and fraction histogram."""
    plt, np = _load_plotting()
    if plt is None or np is None:
        return

    samples = [r[0] for r in report_rows]
    n_hashes = np.array([r[1] for r in report_rows])
    n_in = np.array([r[2] for r in report_rows])
    n_new = np.array([r[4] for r in report_rows])
    pct = np.array([r[3] for r in report_rows])

    # 1) Stacked bar
    fig, ax = plt.subplots(figsize=(max(6, int(len(samples) * 0.6)), 6))
    ax.bar(samples, n_in, label='shared (in DB/other)')
    ax.bar(samples, n_new, bottom=n_in, label='unique (not in DB/other)')
    ax.set_ylabel('Number of unique hashes')
    ax.set_title('Shared vs Unique hashes per sample')
    ax.legend()
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    out1 = f'{plot_prefix}_shared_unique_per_sample.png'
    fig.savefig(out1, dpi=150)
    plt.close(fig)

    # 2) Pairwise Jaccard heatmap – O2: upper-triangle only, then mirror
    S = len(samples)
    jmat = np.zeros((S, S), dtype=float)
    imat = np.zeros((S, S), dtype=int)
    hash_sets = [set(sample_hash_counts[s].keys()) for s in samples]
    for i in range(S):
        jmat[i, i] = 1.0
        imat[i, i] = len(hash_sets[i])
        for j in range(i + 1, S):
            inter = len(hash_sets[i] & hash_sets[j])
            union = len(hash_sets[i] | hash_sets[j])
            val = inter / union if union > 0 else 0.0
            jmat[i, j] = jmat[j, i] = val
            imat[i, j] = imat[j, i] = inter
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.imshow(jmat, cmap='viridis', vmin=0, vmax=1)
    ax.set_xticks(np.arange(S))
    ax.set_yticks(np.arange(S))
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.set_yticklabels(samples)
    ax.set_title('Pairwise Jaccard similarity (hashed k-mers)')
    fig.colorbar(c, ax=ax, label='Jaccard')
    for i in range(S):
        for j in range(S):
            ax.text(j, i, str(imat[i, j]), ha='center', va='center',
                    color='w' if jmat[i, j] < 0.5 else 'black', fontsize=8)
    plt.tight_layout()
    out2 = f'{plot_prefix}_pairwise_jaccard.png'
    fig.savefig(out2, dpi=150)
    plt.close(fig)

    # 3) Histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(pct, bins=min(20, max(4, len(pct))), color='tab:blue', edgecolor='k')
    ax.set_xlabel('Fraction of sample hashes present in DB/other')
    ax.set_ylabel('Number of samples')
    ax.set_title('Distribution of fraction shared')
    plt.tight_layout()
    out3 = f'{plot_prefix}_pct_shared_hist.png'
    fig.savefig(out3, dpi=150)
    plt.close(fig)

    label = 'DB' if db_loaded else 'other samples'
    logger.info('Plots written (compared to %s): %s, %s, %s',
                label, out1, out2, out3)


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser (D2: extracted from main)
# ─────────────────────────────────────────────────────────────────────────────

def build_arg_parser() -> argparse.ArgumentParser:
    """Build and return the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            'Sketch FASTQ, FASTA, or sourmash signature inputs and report '
            'cohort-relative k-mer novelty.'
        )
    )
    parser.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')

    # ── Required ──────────────────────────────────────────────────────────────
    req = parser.add_argument_group('Required parameters')
    req.add_argument('--input-dir', required=True,
                     help='Directory containing input files')
    req.add_argument('--input-type',
                     choices=['paired-fastq', 'single-fastq', 'fasta',
                              'signature'],
                     default='paired-fastq', required=True,
                     help='Type of files in --input-dir (default: paired-fastq)')
    req.add_argument('--out-dir', required=True,
                     help='Directory for all output files')

    # D3: no longer required=True
    req.add_argument('--input-substring', default=None,
                     help='Only filenames containing this substring are considered')

    # B1: --pair was used throughout main but never registered
    req.add_argument('--pair', nargs=3, action='append',
                     metavar=('NAME', 'R1', 'R2'), default=None,
                     help='Explicit sample triple (may be repeated)')

    # ── Sketching ─────────────────────────────────────────────────────────────
    sketch = parser.add_argument_group('Sketching parameters')
    sketch.add_argument('--ksize', type=int, default=31,
                        help='k-mer size (default: 31)')
    sketch.add_argument('--scaled', type=int, default=10000,
                        help='sourmash scaled parameter (default: 10000)')
    sketch.add_argument('--seed', type=int, default=42,
                        help='sourmash hash seed (default: 42)')
    sketch.add_argument('--cache-dir',
                        help='Directory for reusable per-sample sketch caches '
                             '(default: <out-dir>/sketch_cache)')
    sketch.add_argument('--no-cache', action='store_true',
                        help='Disable sketch cache reads/writes')

    # ── Sample metadata ──────────────────────────────────────────────────────
    meta = parser.add_argument_group('Sample metadata')
    meta.add_argument('--metadata',
                      help='Sample metadata CSV/TSV with one row per sample')
    meta.add_argument('--metadata-sample-column', default='sample',
                      help='Column in --metadata containing sample names '
                           '(default: sample)')

    # ── Output files ──────────────────────────────────────────────────────────
    out_grp = parser.add_argument_group('Output file paths')
    # B1: --output was used but never registered
    out_grp.add_argument('--output', default='kmer_counts.csv',
                         help='Path for main CSV output (default: kmer_counts.csv)')
    # B1: --report was used but never registered
    out_grp.add_argument('--report', default='kmer_report.csv',
                         help='Path for per-sample hash-comparison report CSV')
    # B1: --shared-summary was used but never registered
    out_grp.add_argument('--shared-summary', default=None,
                         dest='shared_summary',
                         help='Path for top-shared-per-sample summary CSV')

    # ── Export / limits ───────────────────────────────────────────────────────
    exp = parser.add_argument_group('Export options')
    exp.add_argument('--top-per-sample', type=int, default=200,
                     help='Max top k-mers exported per sample (default: 200)')
    exp.add_argument('--force-top', action='store_true',
                     help='Fail if a sample has fewer than --top-per-sample '
                          'available hashes')
    exp.add_argument('--export-unique-shared', action='store_true',
                     help='Also write per-sample unique/shared hash files')
    exp.add_argument('--fill-unique', action='store_true',
                     help='Pad per-sample unique file to --top-per-sample with '
                          'shared hashes')
    exp.add_argument('--top-shared', type=int, default=100,
                     help='Top shared k-mers per sample (default: 100)')
    exp.add_argument('--long', action='store_true',
                     help='Write long/sparse CSV instead of wide matrix')
    exp.add_argument('--full-exports', action='store_true',
                     help='Enable heavy per-sample hash/sequence exports')

    # ── Profile / DB ──────────────────────────────────────────────────────────
    prof = parser.add_argument_group('Profile and DB parameters')
    prof.add_argument('--profile-load',
                      help='Path to an existing profile directory')
    prof.add_argument('--profile-name',
                      help='Name for this run\'s profile')
    prof.add_argument('--profile-update', action='store_true',
                      help='Update (or create) the named profile')
    prof.add_argument('--profile-force-recompute', action='store_true',
                      help='Recompute samples already in the loaded profile')
    prof.add_argument('--incremental', action='store_true',
                      help='Compare each sample to the DB incrementally')
    prof.add_argument('--db-in',
                      help='Existing hash DB (newline-separated ints, optionally .gz)')
    prof.add_argument('--db-out',
                      help='Path to write hash DB (.gz for gzip)')
    prof.add_argument('--update-db', action='store_true',
                      help='Update the DB with new hashes before writing')
    # B1: --write-db-each was used at line 823 but never registered
    prof.add_argument('--write-db-each', action='store_true',
                      help='Write the DB after processing each sample '
                           '(only with --incremental)')

    # ── K-mer sequence export ─────────────────────────────────────────────────
    kmer = parser.add_argument_group('K-mer sequence export (heavy)')
    kmer.add_argument('--export-kmers-dir',
                      help='Directory for per-sample k-mer CSVs')
    kmer.add_argument('--export-global-shared',
                      help='Filename for global shared k-mers CSV')
    kmer.add_argument('--export-mode', choices=['full', 'top'], default='full',
                      help='Base unique/shared on full MinHash or top-N '
                           '(default: full)')
    kmer.add_argument('--export-sequences', action='store_true', default=False,
                      help='Export canonical k-mer sequences (very heavy)')
    kmer.add_argument('--map-out', default=None,
                      help='Path for hash->sequence mapping CSV')
    kmer.add_argument('--collect-sequences', action='store_true', default=False,
                      help='Collect example k-mer sequences for hashes (slow)')
    kmer.add_argument('--map-threads', type=int, default=4,
                      help='Threads for sequence mapping (default: 4)')
    kmer.add_argument('--global-top', type=int, default=0,
                      help='Limit wide CSV to top-M global hashes (0 = no limit)')
    kmer.add_argument('--export-top-union', action='store_true',
                      help='Write union summary of per-sample top-N hashes')
    kmer.add_argument('--top-union-summary',
                      help='Path for top-union summary CSV')

    # ── Plotting ──────────────────────────────────────────────────────────────
    plot = parser.add_argument_group('Plotting')
    plot.add_argument('--plot', action='store_true',
                      help='Generate shared/unique and Jaccard plots')
    plot.add_argument('--plot-prefix',
                      help='Prefix for plot filenames')

    # ── Outlier detection ─────────────────────────────────────────────────────
    outl = parser.add_argument_group('Outlier detection')
    # B1: all four --outlier-* flags were used but never registered
    outl.add_argument('--export-outlier-only', action='store_true',
                      help='Export hashes unique to outlier samples')
    outl.add_argument('--outlier-threshold', type=float, default=0.5,
                      help='Samples with pct_in_db <= this are outliers '
                           '(default: 0.5)')
    outl.add_argument('--outlier-top', type=int, default=0,
                      help='Limit outlier-only hashes to top-N by count '
                           '(0 = no limit)')
    outl.add_argument('--outlier-dir',
                      help='Directory for outlier output files')
    outl.add_argument('--outlier-seq-per-hash', action='store_true',
                      help='Write example k-mer sequences for outlier hashes')

    # ── Misc ──────────────────────────────────────────────────────────────────
    misc = parser.add_argument_group('Misc')
    misc.add_argument('--verbose', '-v', action='store_true',
                      help='Verbose progress output')
    misc.add_argument('--threads', type=int, default=4,
                      help='Threads for parallel sample sketching (default: 4)')
    misc.add_argument('--allow-partial', action='store_true',
                      help='Continue when one or more samples fail processing')
    misc.add_argument('--multiqc', action='store_true',
                      help='Write MultiQC custom-content JSON')
    misc.add_argument('--ro-crate', action='store_true',
                      help='Write minimal ro-crate-metadata.json provenance')

    return parser


def _validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.ksize <= 0:
        parser.error('--ksize must be > 0')
    if args.scaled <= 0:
        parser.error('--scaled must be > 0')
    if args.seed < 0:
        parser.error('--seed must be >= 0')
    if args.threads <= 0:
        parser.error('--threads must be > 0')
    if args.map_threads <= 0:
        parser.error('--map-threads must be > 0')
    if args.top_per_sample < 0:
        parser.error('--top-per-sample must be >= 0')
    if args.top_shared < 0:
        parser.error('--top-shared must be >= 0')
    if args.global_top < 0:
        parser.error('--global-top must be >= 0')
    if args.outlier_top < 0:
        parser.error('--outlier-top must be >= 0')
    if args.outlier_threshold < 0 or args.outlier_threshold > 1:
        parser.error('--outlier-threshold must be between 0 and 1')
    if args.profile_update and not args.profile_name:
        parser.error('--profile-update requires --profile-name')
    if args.write_db_each and not args.incremental:
        parser.error('--write-db-each requires --incremental')
    if args.no_cache and args.cache_dir:
        parser.error('--no-cache cannot be combined with --cache-dir')
    if args.pair and args.input_type != 'paired-fastq':
        parser.error('--pair is only supported with --input-type paired-fastq')
    if args.fill_unique and not args.export_unique_shared:
        args.export_unique_shared = True
    if args.input_type == 'signature':
        if args.collect_sequences or args.export_sequences or args.outlier_seq_per_hash:
            parser.error(
                'signature input does not contain sequences; disable '
                '--collect-sequences, --export-sequences, and '
                '--outlier-seq-per-hash')


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = build_arg_parser()
    argv_list = list(argv) if argv is not None else None
    args = parser.parse_args(argv_list)
    _validate_args(parser, args)

    if sourmash is None:
        detail = f': {SOURMASH_IMPORT_ERROR}' if SOURMASH_IMPORT_ERROR else ''
        parser.error(
            'sourmash could not be imported'
            f'{detail}. Install a compatible sourmash build.')

    # ── Logging setup ─────────────────────────────────────────────────────────
    # D1: builtins.print monkey-patch removed; all output goes through the logger
    logger = logging.getLogger('OutliMer')
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        fmt = logging.Formatter(
            '%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    # ── Output directory ──────────────────────────────────────────────────────
    out_dir: str = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    # Resolve relative output paths under out_dir
    def _resolve(path: Optional[str], default_name: str) -> str:
        if not path:
            return os.path.join(out_dir, default_name)
        if not os.path.isabs(path) and not os.path.dirname(path):
            return os.path.join(out_dir, path)
        return path

    args.output = _resolve(args.output, 'kmer_counts.csv')
    if args.cache_dir:
        args.cache_dir = _resolve(args.cache_dir, os.path.basename(args.cache_dir))
    elif not args.no_cache:
        args.cache_dir = os.path.join(out_dir, 'sketch_cache')
    if args.report:
        args.report = _resolve(args.report, 'report.csv')
    if args.shared_summary:
        args.shared_summary = _resolve(
            args.shared_summary, os.path.basename(args.shared_summary))
    if args.db_out:
        args.db_out = _resolve(args.db_out, os.path.basename(args.db_out))
    if args.export_kmers_dir:
        args.export_kmers_dir = _resolve(args.export_kmers_dir, 'export_kmers')
    else:
        args.export_kmers_dir = os.path.join(out_dir, 'export_kmers')
    if args.map_out:
        args.map_out = _resolve(args.map_out, os.path.basename(args.map_out))
    if args.top_union_summary:
        args.top_union_summary = _resolve(
            args.top_union_summary, os.path.basename(args.top_union_summary))
    if not args.plot_prefix:
        args.plot_prefix = os.path.join(out_dir, 'kmer_plots')
    elif not os.path.isabs(args.plot_prefix) and not os.path.dirname(args.plot_prefix):
        args.plot_prefix = os.path.join(out_dir, args.plot_prefix)

    # ── Gather samples ────────────────────────────────────────────────────────
    logger.info('Scanning %s (type=%s, substring=%s)',
                args.input_dir, args.input_type, args.input_substring)
    try:
        gathered = gather_samples_from_dir(
            args.input_dir, args.input_type, args.input_substring)
    except Exception as exc:
        parser.error(str(exc))

    # Merge programmatically supplied --pair triples
    pairs: List[Tuple[str, str, str]] = []
    if args.pair:
        for entry in args.pair:
            pairs.append(tuple(entry))  # type: ignore[arg-type]
    for name, r1, r2 in gathered:
        pairs.append((name, r1, r2))

    if not pairs:
        parser.error('No input samples found. Check --input-dir / --pair.')
    try:
        validate_sample_inputs(pairs, args.input_type)
    except ValueError as exc:
        parser.error(str(exc))

    logger.info('GATHERED %d sample(s)', len(pairs))
    for name, r1, r2 in pairs:
        if r2:
            logger.info('  %s: %s + %s', name,
                        os.path.basename(r1), os.path.basename(r2))
        else:
            logger.info('  %s: %s (single)', name, os.path.basename(r1))

    metadata: Dict[str, Dict[str, str]] = {}
    metadata_columns: List[str] = []
    if args.metadata:
        try:
            metadata, metadata_columns = load_sample_metadata(
                args.metadata, args.metadata_sample_column)
        except Exception as exc:
            parser.error(f'Failed to load --metadata {args.metadata}: {exc}')
        pair_names = {name for name, _, _ in pairs}
        missing_metadata = sorted(pair_names - set(metadata))
        if missing_metadata:
            logger.warning('Metadata missing for %d sample(s): %s',
                           len(missing_metadata), ','.join(missing_metadata[:10]))
        extra_metadata = sorted(set(metadata) - pair_names)
        if extra_metadata:
            logger.warning('Metadata has %d sample(s) not in this run',
                           len(extra_metadata))

    # ── Load DB ───────────────────────────────────────────────────────────────
    db_hashes: Set[int] = set()
    db_loaded = False
    if args.db_in:
        logger.info('Loading DB from %s …', args.db_in)
        try:
            db_hashes = load_hash_db(args.db_in)
            db_loaded = True
            logger.info('Loaded %d hashes from DB', len(db_hashes))
        except Exception as exc:
            parser.error(f'Failed to load --db-in {args.db_in}: {exc}')

    # ── Load profile ──────────────────────────────────────────────────────────
    profile_dir: Optional[str] = None
    profile_agg_counts: Dict[int, int] = {}
    profile_samples: Set[str] = set()
    if args.profile_load:
        profile_dir = args.profile_load
        if not os.path.isdir(profile_dir):
            parser.error(f'--profile-load is not a directory: {profile_dir}')
        try:
            db_path = os.path.join(profile_dir, 'db.txt.gz')
            if not os.path.exists(db_path):
                db_path = os.path.join(profile_dir, 'db.txt')
            if os.path.exists(db_path):
                db_hashes = load_hash_db(db_path)
                db_loaded = True
            agg_path = os.path.join(profile_dir, 'agg_counts.csv')
            if os.path.exists(agg_path):
                with open(agg_path, 'r', newline='') as af:
                    for row in csv.reader(af):
                        if not row:
                            continue
                        try:
                            h, c = int(row[0]), int(row[1])
                            profile_agg_counts[h] = (
                                profile_agg_counts.get(h, 0) + c)
                        except (ValueError, IndexError) as exc:
                            logger.debug('Skipping bad agg row: %s', exc)
            samples_path = os.path.join(profile_dir, 'samples.txt')
            if os.path.exists(samples_path):
                with open(samples_path) as sf:
                    for line in sf:
                        n = line.strip()
                        if n:
                            profile_samples.add(n)
            logger.info('Loaded profile: %d samples, %d agg hashes, db=%s',
                        len(profile_samples), len(profile_agg_counts), db_loaded)
        except Exception as exc:
            parser.error(f'Failed to load profile from {profile_dir}: {exc}')

    # ── R4: compute top_n once ────────────────────────────────────────────────
    top_n: Optional[int] = (args.top_per_sample
                             if args.top_per_sample > 0 else None)

    # ── Process samples ───────────────────────────────────────────────────────
    (sample_hash_counts,
     sample_files,
     sourmash_seq_map,
     profile_agg_counts,
     profile_samples,
     failures) = process_samples(
        pairs=pairs,
        ksize=args.ksize,
        scaled=args.scaled,
        seed=args.seed,
        input_type=args.input_type,
        collect_sequences=args.collect_sequences,
        top_per_sample=args.top_per_sample,
        profile_samples=profile_samples,
        profile_dir=profile_dir,
        profile_force_recompute=args.profile_force_recompute,
        profile_agg_counts=profile_agg_counts,
        profile_update=args.profile_update,
        profile_name=args.profile_name,
        threads=args.threads,
        cache_dir=None if args.no_cache else args.cache_dir,
        logger=logger,
    )

    if failures and not args.allow_partial:
        logger.error(
            'Stopping because %d sample(s) failed. Use --allow-partial to '
            'continue explicitly.', len(failures))
        return 1
    if failures:
        logger.warning(
            'Continuing with %d successful sample(s); %d failed sample(s)',
            len(sample_hash_counts), len(failures))
    if not sample_hash_counts:
        logger.error('No samples were successfully processed.')
        return 1

    # ── Always write union summary ────────────────────────────────────────────
    try:
        seq_map_arg = sourmash_seq_map if sourmash_seq_map else None
        union_path = write_union_summary(
            sample_hash_counts, out_dir,
            seq_map=seq_map_arg,
            out_path=args.top_union_summary if args.top_union_summary else None)
        logger.info('Wrote union summary to %s', union_path)
    except Exception as exc:
        logger.error('Failed to write union summary: %s', exc)
        return 1

    # ── Per-sample exports ────────────────────────────────────────────────────
    if args.export_kmers_dir and args.full_exports:
        _run_exports(args, sample_hash_counts, sample_files,
                     sourmash_seq_map, top_n, logger)

    # ── Build report rows ─────────────────────────────────────────────────────
    # B7: removed the always-true `'report_rows' in locals()` guard
    report_rows: List[Tuple[str, int, int, float, int]] = []

    if args.incremental:
        # In incremental mode, accumulate comparisons as samples are processed
        for name, _, _ in pairs:
            if name not in sample_hash_counts:
                continue
            sample_hashes = set(sample_hash_counts[name])
            n_in, n_total, frac = compare_sample_to_db(sample_hashes, db_hashes)
            report_rows.append((name, n_total, n_in, round(frac, 6),
                                 n_total - n_in))
            if args.update_db:
                db_hashes.update(sample_hashes)
                db_loaded = True
            if args.db_out and args.write_db_each:
                save_hash_db(args.db_out, db_hashes)
    else:
        for sample in sorted(sample_hash_counts):
            sample_hashes = set(sample_hash_counts[sample])
            if db_loaded:
                n_in, n_total, frac = compare_sample_to_db(
                    sample_hashes, db_hashes)
            else:
                # Compare against union of all other samples
                other = set()
                for s2, counts in sample_hash_counts.items():
                    if s2 != sample:
                        other.update(counts)
                n_in, n_total, frac = compare_sample_to_db(
                    sample_hashes, other)
            report_rows.append((sample, n_total, n_in, round(frac, 6),
                                 n_total - n_in))
            logger.debug('Report %s: total=%d in_db=%d frac=%.4f new=%d',
                         sample, n_total, n_in, frac, n_total - n_in)

    # ── Write report ──────────────────────────────────────────────────────────
    if args.report:
        with open(args.report, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['sample'] + metadata_columns +
                            ['n_hashes', 'n_in_db',
                             'pct_in_db', 'n_new_hashes'])
            for row in report_rows:
                sample = row[0]
                meta_values = [
                    metadata.get(sample, {}).get(column, '')
                    for column in metadata_columns
                ]
                writer.writerow([sample] + meta_values + list(row[1:]))
        logger.info('Wrote report to %s', args.report)

    multiqc_path: Optional[str] = None
    if args.multiqc:
        multiqc_path = write_multiqc_summary(report_rows, out_dir)
        logger.info('Wrote MultiQC custom content to %s', multiqc_path)

    # ── CSV matrix output ─────────────────────────────────────────────────────
    if args.full_exports:
        if args.long:
            out_long = (args.output if args.output.lower().endswith('.csv')
                        else args.output + '.csv')
            if sourmash_seq_map:
                write_long_csv_with_seq(
                    sample_hash_counts, out_long, sourmash_seq_map)
            else:
                write_long_csv(sample_hash_counts, out_long)
            logger.info('Wrote long CSV to %s', out_long)
        else:
            gtop = args.global_top
            if gtop > 0:
                # R6: reuse total_counts helper
                total_counts = _compute_total_counts(sample_hash_counts)
                top_hashes = [h for h, _ in
                              sorted(total_counts.items(), key=lambda x: -x[1])[:gtop]]
                reduced = {s: {h: counts.get(h, 0) for h in top_hashes}
                           for s, counts in sample_hash_counts.items()}
                logger.info('Writing wide CSV (top %d global hashes) to %s',
                            gtop, args.output)
                if sourmash_seq_map:
                    write_wide_csv_with_seq(reduced, args.output,
                                            sourmash_seq_map)
                else:
                    write_wide_csv(reduced, args.output)
            else:
                logger.info('Writing full wide CSV to %s', args.output)
                if sourmash_seq_map:
                    write_wide_csv_with_seq(sample_hash_counts, args.output,
                                            sourmash_seq_map)
                else:
                    write_wide_csv(sample_hash_counts, args.output)
    else:
        logger.info('Minimal mode: skipping detailed CSV output '
                    '(use --full-exports to enable)')

    # ── Outlier detection ─────────────────────────────────────────────────────
    if args.export_outlier_only:
        _detect_outliers(
            args=args,
            report_rows=report_rows,
            sample_hash_counts=sample_hash_counts,
            sample_files=sample_files,
            sourmash_seq_map=sourmash_seq_map,
            db_hashes=db_hashes,
            db_loaded=db_loaded,
            out_dir=out_dir,
            logger=logger,
        )

    # ── Plots ─────────────────────────────────────────────────────────────────
    if args.plot:
        generate_plots(report_rows, sample_hash_counts,
                       args.plot_prefix, db_loaded, logger)

    # ── Save DB ───────────────────────────────────────────────────────────────
    if args.db_out:
        if db_loaded and not args.update_db:
            out_db = db_hashes
        else:
            out_db = set(db_hashes)
            for counts in sample_hash_counts.values():
                out_db.update(counts)
        save_hash_db(args.db_out, out_db)
        logger.info('Wrote DB (%d hashes) to %s', len(out_db), args.db_out)

    # ── Save profile ──────────────────────────────────────────────────────────
    if args.profile_update and args.profile_name:
        prof_out_dir = os.path.join(out_dir, args.profile_name)
        os.makedirs(prof_out_dir, exist_ok=True)
        save_hash_db(os.path.join(prof_out_dir, 'db.txt.gz'),
                     set(db_hashes) | set().union(*sample_hash_counts.values()))
        agg_path = os.path.join(prof_out_dir, 'agg_counts.csv')
        with open(agg_path, 'w', newline='') as af:
            aw = csv.writer(af)
            for h, c in sorted(profile_agg_counts.items()):
                aw.writerow([h, c])
        samples_path = os.path.join(prof_out_dir, 'samples.txt')
        with open(samples_path, 'w') as sf:
            for n in sorted(profile_samples):
                sf.write(n + '\n')
        logger.info('Saved profile "%s" to %s', args.profile_name, prof_out_dir)

    manifest_path = os.path.join(out_dir, 'run_manifest.json')
    manifest = {
        'schema_version': 2,
        'created_at': datetime.now(timezone.utc).isoformat(),
        'command': ['outlimer', *(argv_list if argv_list is not None
                                  else sys.argv[1:])],
        'software_versions': _runtime_versions(),
        'input_type': args.input_type,
        'ksize': args.ksize,
        'scaled': args.scaled,
        'seed': args.seed,
        'n_samples': len(sample_hash_counts),
        'samples': sorted(sample_hash_counts),
        'failed_samples': failures,
        'input_fingerprints': {
            name: _source_fingerprint((r1, r2))
            for name, r1, r2 in pairs
        },
        'metadata_columns': metadata_columns,
        'cache_dir': None if args.no_cache else args.cache_dir,
        'outputs': {
            'report': args.report,
            'union_summary': union_path if 'union_path' in locals() else None,
            'matrix': args.output if args.full_exports else None,
        },
    }
    with open(manifest_path, 'w') as mf:
        json.dump(manifest, mf, indent=2, sort_keys=True)
    logger.info('Wrote run manifest to %s', manifest_path)

    if args.ro_crate:
        input_paths = {
            'input_dir': args.input_dir,
            'metadata': args.metadata or '',
            'db_in': args.db_in or '',
        }
        output_paths = {
            'report': args.report,
            'union_summary': union_path if 'union_path' in locals() else None,
            'matrix': args.output if args.full_exports else None,
            'manifest': manifest_path,
            'multiqc': multiqc_path,
        }
        crate_path = write_ro_crate_metadata(
            out_dir,
            input_paths,
            output_paths,
            {
                'input_type': args.input_type,
                'ksize': args.ksize,
                'scaled': args.scaled,
                'seed': args.seed,
                'outlimer_version': __version__,
                'n_samples': len(sample_hash_counts),
            },
        )
        logger.info('Wrote RO-Crate metadata to %s', crate_path)

    logger.info('Done.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
