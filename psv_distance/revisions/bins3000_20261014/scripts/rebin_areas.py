#!/usr/bin/env python3
"""Repartition existing >1 kb reference-base areas at 3 kb, preserving all inputs.

Only tail-bearing physical copies are queried. Each copy's full usable PSV set
is searched, including markers outside the confident/lowmap/coverage fragments.
Run production work in Slurm; standard-library-only dependencies are imported
from the original analysis scripts without modifying them.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

SCOPES = ('Segdup', 'Segdup_lowmap')
OLD_TAIL = ('1001-5000', '>5000')
NEW_TAIL = ('1001-3000', '>3000')
KEEP_BINS = ('0', '1-50', '51-150', '151-300', '301-1000', 'no_usable_psv')
NEW_ORDER = KEEP_BINS[:5] + NEW_TAIL + KEEP_BINS[5:]


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical(bed, merge):
    return {chrom: vals for chrom, rows in bed.items() if (vals := merge(rows))}


def tail_split(scoped, psvs, copy_start, copy_end, merge, intersect, subtract):
    """Base distance <=3000 is a half-open PSV interval expanded by 3000."""
    expanded = merge((max(copy_start, s - 3000), min(copy_end, e + 3000))
                     for s, e in psvs)
    near = intersect(scoped, expanded)
    return near, subtract(scoped, near)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True,
                        help='New directory below ANALYSIS_ROOT/revisions; never original prepared/.')
    args = parser.parse_args()
    root = args.analysis_root.resolve()
    prepared = root / 'prepared'
    output = args.output.resolve()
    if root / 'revisions' not in output.parents:
        parser.error('--output must be a new derivative directory under ANALYSIS_ROOT/revisions')
    sys.path.insert(0, str(root / 'scripts'))
    from intervals import (genomic_intersect, genomic_subtract, intersect,
                           merge, read_bed, subtract, total_length, write_bed)
    from prepare_scopes import unique_coverage
    from psv_distance import PSVIndex, span_distance

    # Focused boundary/overlap checks; one SNP PSV at 1-based coordinate 5001.
    test_psv = [(5000, 5001)]
    near, far = tail_split([(0, 10000)], test_psv, 0, 10000,
                          merge, intersect, subtract)
    assert near == [(2000, 8001)] and far == [(0, 2000), (8001, 10000)]
    assert span_distance((2000, 2001), test_psv[0]) == 3000
    assert span_distance((1999, 2000), test_psv[0]) == 3001
    assert span_distance((8000, 8001), test_psv[0]) == 3000
    assert span_distance((8001, 8002), test_psv[0]) == 3001
    # Coverage fragmentation must not restrict a PSV's distance search.
    fragmented, _ = tail_split([(2100, 2200), (7900, 8100)], test_psv,
                               0, 10000, merge, intersect, subtract)
    assert fragmented == [(2100, 2200), (7900, 8001)]

    input_paths = [prepared / 'scope_distance_areas.tsv']
    scopes = {scope: read_bed(prepared / f'{scope}.bed') for scope in SCOPES}
    old_beds = {}
    for scope in SCOPES:
        input_paths.append(prepared / f'{scope}.bed')
        for label in OLD_TAIL:
            path = prepared / f'{scope}_distance_{label.replace(">", "gt")}.bed'
            input_paths.append(path)
            old_beds[scope, label] = read_bed(path)
    before = {str(p.relative_to(root)): {'sha256': digest(p), 'size': p.stat().st_size}
              for p in input_paths}
    with open(prepared / 'scope_distance_areas.tsv', newline='') as handle:
        old_rows = list(csv.DictReader(handle, delimiter='\t'))
    old_areas = {(row['scope'], row['distance_bin']): int(row['evaluable_bp'])
                 for row in old_rows}
    assert len(old_areas) == len(old_rows), 'Duplicate original area keys'
    for scope in SCOPES:
        assert {b for s, b in old_areas if s == scope} == set(KEEP_BINS + OLD_TAIL)
        for label in OLD_TAIL:
            assert total_length(old_beds[scope, label]) == old_areas[scope, label]
        assert sum(n for (s, _), n in old_areas.items() if s == scope) == total_length(scopes[scope])
        assert total_length(genomic_intersect(old_beds[scope, OLD_TAIL[0]],
                                             old_beds[scope, OLD_TAIL[1]])) == 0
    old_tail = {}
    for scope in SCOPES:
        chroms = set(old_beds[scope, OLD_TAIL[0]]) | set(old_beds[scope, OLD_TAIL[1]])
        old_tail[scope] = canonical({c: old_beds[scope, OLD_TAIL[0]].get(c, []) +
                                    old_beds[scope, OLD_TAIL[1]].get(c, []) for c in chroms}, merge)
    assert canonical(genomic_intersect(old_tail['Segdup'], scopes['Segdup_lowmap']), merge) == old_tail['Segdup_lowmap']

    autosomes = {'chr' + str(i) for i in range(1, 23)}
    index = PSVIndex(prepared / 'psv_index.sqlite', autosomes)
    coverage = [tuple(row) for row in index.connection.execute(
        'SELECT chrom,start0,end0,family,copy_id FROM coverage') if row['chrom'] in autosomes]
    unique, ambiguous, bycopy = unique_coverage(coverage)
    assert total_length(genomic_subtract(old_tail['Segdup'], unique)) == 0
    assert total_length(genomic_intersect(old_tail['Segdup'], ambiguous)) == 0
    new_near, new_far, assigned = defaultdict(list), defaultdict(list), defaultdict(list)
    copies_used = []
    for copy_id, covered in bycopy.items():
        copy = index.copies[copy_id]
        chrom, start, end = copy['chrom'], copy['start0'], copy['end0']
        tail = intersect(covered, old_tail['Segdup'].get(chrom, []))
        if not tail:
            continue
        assert all(start <= s < e <= end for s, e in tail)
        psvs = merge(row[:2] for row in index.psvis(copy_id).rows)
        if not psvs:
            raise ValueError('Original usable tail assigned to a copy with no PSV: ' + copy_id)
        near, far = tail_split(tail, psvs, start, end, merge, intersect, subtract)
        new_near[chrom].extend(near)
        new_far[chrom].extend(far)
        assigned[chrom].extend(tail)
        copies_used.append({'copy_id': copy_id, 'family': copy['family'],
                            'tail_bp': sum(e - s for s, e in tail),
                            'n_merged_full_copy_psv_spans': len(psvs)})
        # The index cache is in memory only; release each finished physical copy.
        index._sets.pop(copy_id, None)
    index.connection.close()
    assert canonical(assigned, merge) == old_tail['Segdup']
    assert sum(x['tail_bp'] for x in copies_used) == total_length(old_tail['Segdup']), 'Duplicate copy ownership'

    bands = {NEW_TAIL[0]: canonical(new_near, merge), NEW_TAIL[1]: canonical(new_far, merge)}
    derived_beds, new_rows, scope_checks = {}, [], {}
    for scope in SCOPES:
        for label in NEW_TAIL:
            derived_beds[scope, label] = canonical(genomic_intersect(bands[label], scopes[scope]), merge)
        near, far = (derived_beds[scope, label] for label in NEW_TAIL)
        assert total_length(genomic_intersect(near, far)) == 0
        combined = canonical({c: near.get(c, []) + far.get(c, []) for c in set(near) | set(far)}, merge)
        assert combined == old_tail[scope], 'Tail coordinates or coverage changed'
        assert total_length(genomic_subtract(near, old_beds[scope, OLD_TAIL[0]])) == 0
        assert total_length(genomic_subtract(old_beds[scope, OLD_TAIL[1]], far)) == 0
        for label in NEW_ORDER:
            value = total_length(derived_beds[scope, label]) if label in NEW_TAIL else old_areas[scope, label]
            new_rows.append({'scope': scope, 'distance_bin': label, 'evaluable_bp': value})
        n_scope = total_length(scopes[scope])
        n_tail = total_length(old_tail[scope])
        assert sum(row['evaluable_bp'] for row in new_rows if row['scope'] == scope) == n_scope
        assert total_length(near) + total_length(far) == n_tail
        scope_checks[scope] = {'original_tail_bp': n_tail, '1001-3000_bp': total_length(near),
                               '>3000_bp': total_length(far), 'scope_bp': n_scope,
                               'exact_tail_coordinate_conservation': True,
                               'disjoint_tail_bins': True, 'full_scope_area_conserved': True,
                               'first_five_and_no_usable_psv_areas_unchanged': True}
    assert canonical(genomic_intersect(derived_beds['Segdup', NEW_TAIL[0]], scopes['Segdup_lowmap']), merge) == derived_beds['Segdup_lowmap', NEW_TAIL[0]]
    assert canonical(genomic_intersect(derived_beds['Segdup', NEW_TAIL[1]], scopes['Segdup_lowmap']), merge) == derived_beds['Segdup_lowmap', NEW_TAIL[1]]

    # No outputs are written until all geometry and conservation assertions pass.
    output.mkdir(parents=True, exist_ok=True)
    for (scope, label), bed in derived_beds.items():
        write_bed(output / f'{scope}_distance_{label.replace(">", "gt")}.bed', bed)
    with open(output / 'scope_distance_areas.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['scope', 'distance_bin', 'evaluable_bp'], delimiter='\t')
        writer.writeheader()
        writer.writerows(new_rows)
    after = {str(p.relative_to(root)): {'sha256': digest(p), 'size': p.stat().st_size}
             for p in input_paths}
    assert before == after, 'Original prepared inputs changed during execution'
    report = {'validated': True, 'timestamp_utc': datetime.now(timezone.utc).isoformat(),
              'analysis_root': str(root), 'output': str(output), 'original_inputs_unchanged': True,
              'original_inputs': before, 'new_bins': list(NEW_ORDER),
              'scopes': scope_checks, 'tail_copies_queried': len(copies_used),
              'tail_families': len({x['family'] for x in copies_used}),
              'copies_queried': copies_used, 'focused_boundary_checks_passed': True,
              'definition': 'Reference-base distance to usable actual PSV mismatch runs or trimmed span/anchor in the same physical copy. Full-copy PSV search; only the old >1000 bp base-area union is repartitioned.',
              'limitation': 'Reference-base area is not an indel-span exposure denominator.',
              'imports_sha256': {name: digest(root / 'scripts' / name) for name in
                                  ('intervals.py', 'psv_distance.py', 'prepare_scopes.py', 'prepare_psv_database.py')}}
    (output / 'area_rebin_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'validated': True, 'scopes': scope_checks,
                      'tail_copies_queried': len(copies_used)}, indent=2))


if __name__ == '__main__':
    main()
