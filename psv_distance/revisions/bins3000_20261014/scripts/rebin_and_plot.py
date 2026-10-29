#!/usr/bin/env python3
"""Re-bin existing official observations; preserve matching, filtering and distances."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import zipfile
import numpy as np
import pandas as pd
import summary_core_original as core
from plot_clean import BINS, plot

OLD_BINS = list(core.MAIN_BINS)
OLD_ALL = list(core.BIN_ORDER)
GROUP = core.GROUP
KEY = GROUP + ['vtype', 'distance_bin']
COUNTS = core.COUNTS


def distance_bin(distance, cutoff):
    if distance < 0:
        raise ValueError('Negative distance')
    for limit, label in [(0,'0'), (50,'1-50'), (150,'51-150'),
                         (300,'151-300'), (1000,'301-1000'),
                         (cutoff,f'1001-{cutoff}')]:
        if distance <= limit:
            return label
    return f'>{cutoff}'


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20260929)
    a = parser.parse_args()
    root, out = a.analysis_root.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    assert a.bootstrap > 0
    assert [distance_bin(d,3000) for d in [0,1,50,51,150,151,300,301,1000,1001,3000,3001,5000,5001]] == [
        '0','1-50','1-50','51-150','51-150','151-300','151-300','301-1000','301-1000',
        '1001-3000','1001-3000','>3000','>3000','>3000']
    paths = [root/'results'/f'{tech}_DeepVariant_{vtype}.sites.tsv.gz'
             for tech in ['SR','LR'] for vtype in ['SNP','INDEL']]
    source_metrics = root/'results/summary/metrics_DeepVariant.tsv'
    all_sources = paths + [source_metrics, root/'prepared/scope_distance_areas.tsv',
                          root/'scripts/summarize_plot.py']
    before = {str(p):file_hash(p) for p in all_sources}
    family_counts = defaultdict(Counter)
    old_counts = defaultdict(Counter)
    seen = set()
    migrated = Counter()
    count_values = Counter()
    rows_read = 0
    assignment_fields = ['source_file','callset','tech','caller','vtype','scope','observation_id',
                         'side','all_decision','pass_decision','count','distance_bp',
                         'state','family','old_distance_bin','distance_bin']
    with gzip.open(out/'site_bin_assignments.tsv.gz', 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=assignment_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for path in paths:
            with gzip.open(path,'rt') as source:
                for row in csv.DictReader(source,delimiter='\t'):
                    ident = (row['callset'],row['vtype'],row['scope'],row['observation_id'])
                    if ident in seen:
                        raise ValueError('Duplicate side observation: '+repr(ident))
                    seen.add(ident)
                    assert row['caller']=='DeepVariant'
                    assert row['side'] in {'TRUTH','QUERY'}
                    count = int(row['count'])
                    assert count == 1, 'Unexpected official observation weight; investigate before rebinning.'
                    count_values[count] += 1
                    old = row['distance_bin']
                    if row['state']=='usable':
                        d = int(row['distance_bp'])
                        assert old == distance_bin(d,5000)
                        new = distance_bin(d,3000)
                    else:
                        assert old in core.SPECIAL_BINS
                        new = old
                    migrated[old,new] += count
                    for policy, field in [('ALL','all_decision'),('PASS','pass_decision')]:
                        decision = row[field]
                        klass = row['side']+'.'+decision
                        if decision == 'EXCLUDED_FILTER':
                            assert row['side']=='QUERY'
                            klass = 'QUERY.FILTERED'
                        if klass not in COUNTS:
                            raise ValueError('Unexpected benchmark class: '+klass)
                        prefix = (row['callset'],row['tech'],row['caller'],row['scope'],policy,row['vtype'])
                        family_counts[prefix+(new,row['family'])][klass] += count
                        old_counts[prefix+(old,)][klass] += count
                    annotation = {key:row[key] for key in assignment_fields if key in row}
                    annotation.update(source_file=str(path.relative_to(root)),old_distance_bin=old,distance_bin=new)
                    writer.writerow(annotation)
                    rows_read += 1
    old_frame = pd.read_csv(source_metrics,sep='\t',dtype={'distance_bin':str})
    assert not old_frame.duplicated(KEY).any()
    old_expected = {tuple(row[k] for k in KEY): {c:int(row[c]) for c in COUNTS}
                    for row in old_frame.to_dict('records')}
    assert not (set(old_counts)-set(old_expected))
    for key, counts in old_expected.items():
        assert all(old_counts[key][c]==counts[c] for c in COUNTS), key
    print(f'Reconstructed all {len(old_frame)} original metric rows exactly from {rows_read} unique side/scope observations.',flush=True)
    records = [{**dict(zip(KEY+['family'],key)),**{c:value[c] for c in COUNTS}}
               for key,value in sorted(family_counts.items())]
    family_path = out/'family_counts_DeepVariant.tsv'
    pd.DataFrame(records).to_csv(family_path,sep='\t',index=False)
    # Isolated runtime configuration of the exact original bootstrap implementation.
    # The source snapshot and original project scripts are not edited.
    core.MAIN_BINS = list(BINS)
    core.BIN_ORDER = list(BINS)+core.SPECIAL_BINS
    data = core.load_counts([family_path])
    metrics, contrasts = core.summarize(data,a.bootstrap,a.seed)
    metrics = core.add_areas(metrics,out/'prepared/scope_distance_areas.tsv')
    total_keys = GROUP+['vtype']
    old_totals = old_frame.groupby(total_keys)[COUNTS].sum().sort_index()
    new_totals = metrics.groupby(total_keys)[COUNTS].sum().sort_index()
    pd.testing.assert_frame_equal(old_totals,new_totals,check_dtype=False)
    unchanged_bins = BINS[:5] + core.SPECIAL_BINS
    stable_columns = COUNTS+['n_truth','n_query','n_families_truth','n_families_query',
                            'recall','precision','F1','recall_ci_low','recall_ci_high',
                            'precision_ci_low','precision_ci_high','F1_ci_low','F1_ci_high']
    old_stable = old_frame[old_frame.distance_bin.isin(unchanged_bins)].set_index(KEY)[stable_columns].sort_index()
    new_stable = metrics[metrics.distance_bin.isin(unchanged_bins)].set_index(KEY)[stable_columns].sort_index()
    pd.testing.assert_frame_equal(old_stable,new_stable,check_dtype=False,check_exact=False,rtol=1e-8,atol=1e-9)
    for name, frame in [('metrics_DeepVariant.tsv',metrics),
                        ('metrics_PASS_DeepVariant.tsv',metrics[metrics.policy=='PASS'])]:
        frame.to_csv(out/name,sep='\t',index=False,na_rep='NaN',float_format='%.12g')
    # Full support and CIs stay numeric; the eight PNGs contain only the main curves.
    support = metrics[(metrics.policy=='PASS') & metrics.distance_bin.isin(BINS)]
    support.to_csv(out/'support_DeepVariant.tsv',sep='\t',index=False,na_rep='NaN',float_format='%.12g')
    concise = support[['scope','vtype','tech','distance_bin','n_truth','n_query',
                       'n_families_truth','n_families_query']].copy()
    for metric in ['recall','precision','F1']:
        concise[metric+'_percent'] = support[metric]*100
    concise.to_csv(out/'recall_precision_percent.tsv',sep='\t',index=False,na_rep='NaN',float_format='%.5f')
    figures = plot(metrics,out/'figures')
    assert len(figures)==8
    assert not list(out.rglob('*.pdf'))
    zip_path = out/'DeepVariant_PSV_3000bp_PNG.zip'
    with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for figure in figures:
            z.write(figure,Path(figure).name)
    with zipfile.ZipFile(zip_path) as z:
        assert len(z.namelist())==8 and all(p.endswith('.png') for p in z.namelist())
        assert z.testzip() is None
    after = {str(p):file_hash(p) for p in all_sources}
    assert before == after, 'Original input modified during revision'
    import matplotlib
    provenance = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        command=sys.argv,slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__),
        n_observations=rows_read,observation_count_values=dict(count_values),
        old_bins=OLD_BINS,new_bins=BINS,bin_selection='Post hoc user-requested 3 kb revision; original 5 kb results preserved.',
        changed_bin_migrations=[dict(old=x,new=y,n_observations=n) for (x,y),n in sorted(migrated.items()) if x!=y],
        recomputed=['bin assignment','family counts','metrics','95% family-bootstrap intervals','tail reference-base areas','PNG figures'],
        reused=['full-context hap.py matching','official side counts and PASS decisions','normalized distances','PSV catalogue','scopes'],
        original_metrics_reconstructed_exactly=True,overall_counts_conserved=True,
        unchanged_first_five_bins_and_special_states=True,unchanged_bin_intervals_match=True,
        input_hashes_unchanged=True,input_sha256=before,
        bootstrap_replicates=a.bootstrap,bootstrap_seed=a.seed,
        ci_method='Original joint duplication-family percentile bootstrap, 95%; fewer than two contributing families -> undefined CI.',
        sensitivity_note='Other prior sensitivity definitions retain original 5 kb bins; not rerun in this presentation revision.',
        plot_note='No CI or support panel. All defined estimates shown, including n<50; undefined rates remain NaN. Identical categorical x positions, no physical-distance interpolation.',
        figures=figures,png_zip=str(zip_path))
    (out/'validation_and_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps({k:provenance[k] for k in ['n_observations','original_metrics_reconstructed_exactly','overall_counts_conserved','input_hashes_unchanged','png_zip']},indent=2),flush=True)


if __name__=='__main__':
    main()
