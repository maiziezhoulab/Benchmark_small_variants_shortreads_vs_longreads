#!/usr/bin/env python3
"""Pool PSV-distance benchmark counts, bootstrap duplication families, and plot.

Run on Slurm for the real dataset. Counts and confidence intervals are calculated
separately for truth/query sides. PASS is used for figures; ALL is retained as a
numeric sensitivity analysis. The seven distance bins are fixed before analysis.
"""
import argparse
from datetime import datetime, timezone
import glob
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

MAIN_BINS = ['0', '1-50', '51-150', '151-300', '301-1000', '1001-5000', '>5000']
SPECIAL_BINS = ['no_usable_psv', 'ambiguous_copy', 'database_uncovered',
                'representation_crosses_copy', 'unsupported_allele']
BIN_ORDER = MAIN_BINS + SPECIAL_BINS
COUNTS = ['TRUTH.TP', 'TRUTH.FN', 'QUERY.TP', 'QUERY.FP', 'QUERY.UNK', 'QUERY.FILTERED']
GROUP = ['callset', 'tech', 'caller', 'scope', 'policy']
CELL = ['vtype', 'distance_bin']
COLORS = {'DeepVariant': '#0072B2', 'GATK': '#D55E00', 'Strelka2': '#009E73',
          'Clair3': '#CC79A7', 'Longshot': '#E69F00', 'NanoCaller': '#6A3D9A'}
BAD_FAMILIES = {'', '.', 'NA', 'NaN', 'nan', 'None', 'none', 'unassigned',
                'unknown', 'database_uncovered', 'ambiguous_copy', 'unsupported_allele'}


def divide(numerator, denominator):
    numerator, denominator = np.broadcast_arrays(numerator, denominator)
    return np.divide(numerator, denominator,
                     out=np.full(numerator.shape, np.nan, dtype=float),
                     where=denominator != 0)


def count_metrics(counts):
    """The last axis is TPtruth,FNtruth,TPquery,FPquery,UNKquery,filteredquery; the last two never enter denominators."""
    recall = divide(counts[..., 0], counts[..., 0] + counts[..., 1])
    precision = divide(counts[..., 2], counts[..., 2] + counts[..., 3])
    f1 = divide(2 * recall * precision, recall + precision)
    # If both denominators exist and both ratios are zero, F1 is zero.
    f1 = np.where((recall == 0) & (precision == 0), 0., f1)
    return np.stack([recall, precision, f1], axis=-1)


def finite_interval(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return np.nan, np.nan, 0
    low, high = np.quantile(finite, [.025, .975])
    return float(low), float(high), len(finite)


def key_seed(seed, key):
    digest = hashlib.sha256(json.dumps(list(key)).encode()).digest()
    return np.random.SeedSequence([seed, int.from_bytes(digest[:4], 'big')])


def pool_group(frame, key, n_bootstrap, seed):
    families = sorted(frame['family'].unique())
    vtypes = [v for v in ['SNP', 'INDEL'] if v in set(frame['vtype'])]
    cells = [(vtype, distance) for vtype in vtypes for distance in BIN_ORDER]
    family_lookup = {family: i for i, family in enumerate(families)}
    cell_lookup = {cell: i for i, cell in enumerate(cells)}
    cube = np.zeros((len(families), len(cells), len(COUNTS)), dtype=np.int64)
    unknown_family = np.zeros(len(families), dtype=bool)
    for family, i in family_lookup.items():
        unknown_family[i] = family in BAD_FAMILIES or family.startswith(('unassigned:', '__'))
    for row in frame.itertuples(index=False, name=None):
        # Input columns are explicitly arranged by the caller.
        family, vtype, distance, *counts = row
        cube[family_lookup[family], cell_lookup[(vtype, distance)], :] += np.array(counts, dtype=np.int64)
    totals = cube.sum(axis=0)
    points = count_metrics(totals)
    denom_truth = cube[:, :, 0] + cube[:, :, 1]
    denom_query = cube[:, :, 2] + cube[:, :, 3]
    n_truth_families = (denom_truth > 0).sum(axis=0)
    n_query_families = (denom_query > 0).sum(axis=0)
    family_valid_truth = ((denom_truth > 0) & unknown_family[:, None]).sum(axis=0) == 0
    family_valid_query = ((denom_query > 0) & unknown_family[:, None]).sum(axis=0) == 0
    ci_eligible = np.stack([
        (n_truth_families >= 2) & family_valid_truth,
        (n_query_families >= 2) & family_valid_query,
        (n_truth_families >= 2) & (n_query_families >= 2) & family_valid_truth & family_valid_query,
    ], axis=-1)
    draws = np.full((n_bootstrap, len(cells), 3), np.nan)
    # Unassigned-status labels are bookkeeping categories, never biological clusters.
    bootstrap_cube = cube[~unknown_family]
    n_bootstrap_families = len(bootstrap_cube)
    if n_bootstrap_families >= 2 and n_bootstrap > 0:
        rng = np.random.default_rng(key_seed(seed, key))
        flattened = bootstrap_cube.reshape(n_bootstrap_families, -1).astype(float)
        for start in range(0, n_bootstrap, 128):
            batch = min(128, n_bootstrap - start)
            weights = rng.multinomial(n_bootstrap_families, np.full(n_bootstrap_families, 1 / n_bootstrap_families), size=batch)
            sample_counts = (weights @ flattened).reshape(batch, len(cells), len(COUNTS))
            draws[start:start + batch] = count_metrics(sample_counts)
    rows = []
    for cell_index, (vtype, distance) in enumerate(cells):
        nt = int(totals[cell_index, 0] + totals[cell_index, 1])
        nq = int(totals[cell_index, 2] + totals[cell_index, 3])
        row = dict(zip(GROUP, key))
        row.update(vtype=vtype, distance_bin=distance)
        row.update({label: int(totals[cell_index, j]) for j, label in enumerate(COUNTS)})
        row.update(n_truth=nt, n_query=nq, n_query_filtered=int(totals[cell_index, 5]),
                   n_families_truth=int(n_truth_families[cell_index]),
                   n_families_query=int(n_query_families[cell_index]),
                   n_families_union=int(((denom_truth[:, cell_index] + denom_query[:, cell_index]) > 0).sum()),
                   largest_family_fraction_truth=float(denom_truth[:, cell_index].max() / nt) if nt else np.nan,
                   largest_family_fraction_query=float(denom_query[:, cell_index].max() / nq) if nq else np.nan,
                   bootstrap_population_families=n_bootstrap_families)
        for metric_index, metric in enumerate(['recall', 'precision', 'F1']):
            row[metric] = float(points[cell_index, metric_index])
            if ci_eligible[cell_index, metric_index] and n_bootstrap:
                low, high, valid = finite_interval(draws[:, cell_index, metric_index])
            else:
                low, high, valid = np.nan, np.nan, 0
            row[metric + '_ci_low'] = low
            row[metric + '_ci_high'] = high
            row[metric + '_bootstrap_valid'] = valid
        for side, count, nfamily, valid_family in [
                ('truth', nt, n_truth_families[cell_index], family_valid_truth[cell_index]),
                ('query', nq, n_query_families[cell_index], family_valid_query[cell_index])]:
            status = []
            if not count:
                status.append('no_support')
            else:
                if count < 50:
                    status.append('n_lt_50')
                if not valid_family:
                    status.append('family_unassigned')
                if nfamily < 2:
                    status.append('single_family_no_ci')
                elif nfamily < 20:
                    status.append('families_lt_20')
            row['support_status_' + side] = ';'.join(status) or 'adequate_descriptive_support'
        rows.append(row)

    # Descriptive contrasts retain the same joint family draws as the bin CIs.
    comparisons = []
    for vtype in vtypes:
        for metric_index, metric in enumerate(['recall', 'precision']):
            finite_bins = [d for d in MAIN_BINS[1:]
                           if np.isfinite(points[cell_lookup[(vtype, d)], metric_index])]
            if len(finite_bins) < 2:
                continue
            for first, last, contrast_type in [
                    (finite_bins[0], finite_bins[-1], 'farthest_minus_nearest_nonzero'),
                    *[(first, last, 'adjacent_supported_nonzero')
                      for first, last in zip(finite_bins, finite_bins[1:])]]:
                left, right = cell_lookup[(vtype, first)], cell_lookup[(vtype, last)]
                delta = points[right, metric_index] - points[left, metric_index]
                if ci_eligible[left, metric_index] and ci_eligible[right, metric_index]:
                    lo, hi, nfinite = finite_interval(draws[:, right, metric_index] - draws[:, left, metric_index])
                else:
                    lo, hi, nfinite = np.nan, np.nan, 0
                comparisons.append({**dict(zip(GROUP, key)), 'vtype': vtype, 'metric': metric,
                                    'contrast': contrast_type, 'from_bin': first, 'to_bin': last,
                                    'delta': float(delta), 'delta_ci_low': lo,
                                    'delta_ci_high': hi, 'bootstrap_valid': nfinite,
                                    'interpretation': 'Descriptive; no monotonic model or multiplicity-adjusted test.'})
    return rows, comparisons


def load_counts(paths):
    frames = []
    expected = GROUP + CELL + ['family'] + COUNTS
    for path in paths:
        frame = pd.read_csv(path, sep='\t', dtype={column: str for column in GROUP + CELL + ['family']},
                            keep_default_na=False)
        if 'QUERY.FILTERED' not in frame.columns:
            frame['QUERY.FILTERED'] = 0
        missing = set(expected) - set(frame.columns)
        if missing:
            raise ValueError(str(path) + ': missing columns ' + ','.join(sorted(missing)))
        for count in COUNTS:
            frame[count] = pd.to_numeric(frame[count], errors='raise')
            if (frame[count] < 0).any() or ((frame[count] % 1) != 0).any():
                raise ValueError(str(path) + ': noninteger or negative ' + count)
            frame[count] = frame[count].astype(np.int64)
        frames.append(frame[expected])
    if not frames:
        raise ValueError('No family-count files matched.')
    data = pd.concat(frames, ignore_index=True)
    if data.duplicated(GROUP + CELL + ['family']).any():
        raise ValueError('Duplicate callset/scope/policy/type/bin/family keys across input files.')
    for column, allowed in [('scope', {'Segdup', 'Segdup_lowmap'}),
                            ('tech', {'SR', 'LR'}), ('policy', {'PASS', 'ALL'}),
                            ('vtype', {'SNP', 'INDEL'}), ('distance_bin', set(BIN_ORDER))]:
        unknown = set(data[column]) - allowed
        if unknown:
            raise ValueError('Unexpected ' + column + ': ' + repr(unknown))
    return data


def summarize(data, n_bootstrap, seed):
    rows, comparisons = [], []
    for key, frame in data.groupby(GROUP, sort=True, dropna=False):
        frame = frame[['family', 'vtype', 'distance_bin'] + COUNTS]
        local_rows, local_comparisons = pool_group(frame, key, n_bootstrap, seed)
        rows.extend(local_rows)
        comparisons.extend(local_comparisons)
    return pd.DataFrame(rows), pd.DataFrame(comparisons)


def add_areas(metrics, area_path):
    area = pd.read_csv(area_path, sep='\t', dtype={'distance_bin': str})
    required = ['scope', 'distance_bin', 'evaluable_bp']
    if not set(required).issubset(area.columns):
        raise ValueError('Area table must have scope, distance_bin, evaluable_bp.')
    if area.duplicated(['scope', 'distance_bin']).any():
        raise ValueError('Area table has duplicate scope/bin keys.')
    if (pd.to_numeric(area['evaluable_bp']) < 0).any():
        raise ValueError('Negative evaluable area.')
    metrics = metrics.merge(area[required], on=['scope', 'distance_bin'], how='left', validate='many_to_one')
    metrics['area_interpretation'] = 'Single-reference-base geometry; not indel-span exposure.'
    return metrics


def plot_metrics(metrics, output_dir):
    """Present DeepVariant SR and HiFi together; never recalculate estimates or CIs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = []
    labels = ['0\n(overlap)', '1–50', '51–150', '151–300', '301–1,000', '1,001–5,000', '>5,000']
    source = metrics[(metrics['caller'] == 'DeepVariant') &
                     (metrics['policy'] == 'PASS') & metrics['distance_bin'].isin(MAIN_BINS)]
    styles = {'SR': dict(color='#0072B2', marker='o', ls='-', label='SR'),
              'LR': dict(color='#D55E00', marker='s', ls='--', label='HiFi (LR)')}
    for scope in ['Segdup', 'Segdup_lowmap']:
        for vtype in ['SNP', 'INDEL']:
            panel = source[(source['scope'] == scope) & (source['vtype'] == vtype)]
            if panel.empty:
                continue
            if panel.duplicated(['tech', 'distance_bin']).any():
                raise ValueError('More than one DeepVariant observation per technology/bin.')
            if set(panel['tech']) != {'SR', 'LR'}:
                raise ValueError('Both DeepVariant technologies are required in each figure.')
            for metric in ['recall', 'precision']:
                figure, axis = plt.subplots(figsize=(8.8, 5.35))
                for tech in ['SR', 'LR']:
                    values = panel[panel['tech'] == tech].set_index('distance_bin').reindex(MAIN_BINS)
                    x = np.arange(len(MAIN_BINS))
                    y = values[metric].to_numpy(dtype=float) * 100
                    lower = values[metric + '_ci_low'].to_numpy(dtype=float) * 100
                    upper = values[metric + '_ci_high'].to_numpy(dtype=float) * 100
                    style = styles[tech]
                    face = style['color'] if tech == 'SR' else 'white'
                    size = 5.5 if tech == 'SR' else 6.5
                    # Same categorical x in both technologies. NaNs retain line breaks.
                    axis.plot(x[1:], y[1:], color=style['color'], marker=style['marker'],
                              ls=style['ls'], lw=2, ms=size, mfc=face, mew=1.3,
                              label=style['label'], zorder=4 if tech == 'SR' else 3)
                    axis.plot(x[:1], y[:1], color=style['color'], marker=style['marker'],
                              ls='none', ms=size, mfc=face, mew=1.3, zorder=4 if tech == 'SR' else 3)
                    good = np.isfinite(y) & np.isfinite(lower) & np.isfinite(upper)
                    # Direct whiskers allow percentile intervals outside the point estimate.
                    axis.vlines(x[good], lower[good], upper[good], color=style['color'], lw=1.3, alpha=.65, zorder=2)
                    axis.hlines(lower[good], x[good]-.04, x[good]+.04, color=style['color'], lw=1.3, alpha=.65)
                    axis.hlines(upper[good], x[good]-.04, x[good]+.04, color=style['color'], lw=1.3, alpha=.65)
                    if np.any((y[np.isfinite(y)] < 0) | (y[np.isfinite(y)] > 100)):
                        raise ValueError('Metric outside the probability range.')
                axis.set_ylim(-3, 103)
                axis.set_yticks(np.arange(0, 101, 20))
                axis.set_xlim(-.4, len(MAIN_BINS)-.6)
                axis.set_xticks(np.arange(len(MAIN_BINS)), labels, fontsize=10)
                axis.axvline(.5, color='.7', ls=':', lw=1)
                axis.grid(axis='y', alpha=.2)
                axis.set_xlabel('Distance to nearest reference PSV (bp; categorical bins)', labelpad=10)
                axis.set_ylabel(metric.capitalize()+' (%)')
                scope_label = 'Segmental duplications' if scope == 'Segdup' else 'Segmental duplications ∩ low-mappability'
                axis.set_title(f'DeepVariant | {vtype} {metric}\n{scope_label}', fontsize=14, pad=14)
                axis.legend(loc='lower left', frameon=False, fontsize=11)
                figure.tight_layout(pad=1.2)
                stem = output_dir / f'{scope}_{vtype}_{metric}_vs_psv_distance'
                for suffix in ['png', 'pdf']:
                    path = stem.with_suffix('.'+suffix)
                    figure.savefig(path, dpi=220, bbox_inches='tight')
                    files.append(str(path))
                plt.close(figure)
    return files


def write_deepvariant_tables(output_dir):
    """Filter saved source rows verbatim; support details remain outside the figures."""
    import csv
    output_dir = Path(output_dir)
    written = []
    for name in ['metrics', 'overall_scope_metrics', 'sensitivity_comparison']:
        source = output_dir / (name+'.tsv')
        if not source.exists():
            continue
        with source.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            columns = reader.fieldnames
            rows = [r for r in reader if r['caller'] == 'DeepVariant']
        target = output_dir / (name+'_DeepVariant.tsv')
        with target.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
        written.append(str(target))
        if name == 'metrics':
            fields = ['scope', 'vtype', 'tech', 'callset', 'policy', 'distance_bin',
                      'TRUTH.TP', 'TRUTH.FN', 'QUERY.TP', 'QUERY.FP', 'n_truth', 'n_query',
                      'n_families_truth', 'n_families_query', 'evaluable_bp',
                      'recall', 'recall_ci_low', 'recall_ci_high',
                      'precision', 'precision_ci_low', 'precision_ci_high',
                      'support_status_truth', 'support_status_query']
            fields = [c for c in fields if c in columns]
            target = output_dir/'support_DeepVariant.tsv'
            with target.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(r for r in rows if r['policy']=='PASS' and r['distance_bin'] in MAIN_BINS)
            written.append(str(target))
    return written



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family-counts', nargs='+', required=True,
                        help='TSV paths or quoted glob patterns; duplicate paths are removed.')
    parser.add_argument('--areas', help='Primary reference-base area table; omit for sensitivity definitions without recomputed areas.')
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20260929)
    parser.add_argument('--skip-plots', action='store_true')
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error('--bootstrap cannot be negative')
    paths = sorted({str(Path(path).resolve()) for pattern in args.family_counts
                    for path in glob.glob(pattern)})
    data = load_counts(paths)
    metrics, comparisons = summarize(data, args.bootstrap, args.seed)
    if args.areas:
        metrics = add_areas(metrics, args.areas)
    else:
        metrics['evaluable_bp'] = np.nan
        metrics['area_interpretation'] = 'Not recomputed for this sensitivity definition.'
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(output / 'metrics.tsv', sep='\t', index=False, na_rep='NaN', float_format='%.10g')
    metrics[metrics['policy'] == 'PASS'].to_csv(output / 'metrics_PASS.tsv', sep='\t', index=False, na_rep='NaN', float_format='%.10g')
    metrics[metrics['policy'] == 'ALL'].to_csv(output / 'metrics_ALL.tsv', sep='\t', index=False, na_rep='NaN', float_format='%.10g')
    comparisons.to_csv(output / 'descriptive_contrasts.tsv', sep='\t', index=False, na_rep='NaN', float_format='%.10g')
    presentation_tables = write_deepvariant_tables(output)
    figures = [] if args.skip_plots else plot_metrics(metrics, output / 'figures')
    provenance = {'timestamp_utc': datetime.now(timezone.utc).isoformat(), 'bootstrap': args.bootstrap,
                  'seed': args.seed, 'cluster': 'DuploMap connected duplication family',
                  'resampling': 'All observed families jointly across bins/types within callset/scope/policy; percentile 95% CI.',
                  'software': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__},
                  'inputs': [{'path': path, 'bytes': Path(path).stat().st_size,
                              'sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest()} for path in paths],
                  'area_file': str(Path(args.areas).resolve()) if args.areas else None, 'figures': figures,
                  'presentation_caller': 'DeepVariant', 'presentation_tables': presentation_tables,
                  'total_counts_conserved': {column: int(data[column].sum()) == int(metrics[column].sum()) for column in COUNTS},
                  'main_bins': MAIN_BINS, 'special_categories': SPECIAL_BINS}
    if not all(provenance['total_counts_conserved'].values()):
        raise RuntimeError('Grouped count conservation failed.')
    (output / 'summary_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'metric_rows': len(metrics), 'figures': len(figures), 'output': str(output),
                      'conserved': provenance['total_counts_conserved']}))


if __name__ == '__main__':
    main()
