"""PNG-only PSV presentation; categorical positions, no CI or support panel."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

BINS = ['0', '1-50', '51-150', '151-300', '301-1000', '1001-3000', '>3000']
LABELS = ['0\n(overlap)', '1–50', '51–150', '151–300', '301–1,000', '1,001–3,000', '>3,000']


def plot(metrics, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.titleweight': 'normal'})
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source = metrics[(metrics.caller == 'DeepVariant') & (metrics.policy == 'PASS')]
    styles = {'SR': dict(color='#C0392B', marker='o', label='SR'),
              'LR': dict(color='#2463A5', marker='s', label='HiFi (LR)')}
    plotted = []
    files = []
    for scope in ['Segdup', 'Segdup_lowmap']:
        for vtype in ['SNP', 'INDEL']:
            data = source[(source.scope == scope) & (source.vtype == vtype)]
            if data.duplicated(['tech', 'distance_bin']).any():
                raise ValueError('Duplicate technology/bin row.')
            assert set(data.tech) == {'SR', 'LR'}
            for metric in ['recall', 'precision']:
                fig, ax = plt.subplots(figsize=(10.8, 6.3))
                x = np.arange(len(BINS))
                for tech in ['SR', 'LR']:
                    frame = data[data.tech == tech].set_index('distance_bin').reindex(BINS)
                    y = frame[metric].to_numpy(float) * 100
                    finite = y[np.isfinite(y)]
                    assert np.all((finite >= 0) & (finite <= 100))
                    style = styles[tech]
                    # No jitter, midpoint, empirical median, interpolation, or fitted model.
                    kw = dict(color=style['color'], marker=style['marker'], lw=2.2,
                              ms=7 if tech == 'SR' else 7.5, mew=1.5,
                              mfc=style['color'] if tech == 'SR' else 'white',
                              zorder=4 if tech == 'SR' else 3)
                    ax.plot(x[1:], y[1:], ls='-', label=style['label'], **kw)
                    ax.plot(x[:1], y[:1], ls='none', **kw)
                    for i, distance in enumerate(BINS):
                        plotted.append(dict(scope=scope, vtype=vtype, metric=metric,
                                            tech=tech, distance_bin=distance, x=int(x[i]),
                                            value_percent=float(y[i]),
                                            n_truth=int(frame.iloc[i].n_truth),
                                            n_query=int(frame.iloc[i].n_query)))
                # Tiny margins expose markers at 0 and 100; labeled scale is 0–100%.
                ax.set_ylim(-2.5, 102.5)
                ax.set_yticks(np.arange(0, 101, 20))
                ax.set_xlim(-.38, 6.38)
                ax.set_xticks(x, LABELS, fontsize=12)
                ax.tick_params(axis='both', length=4)
                ax.grid(axis='y', color='#d9d9d9', lw=.8, alpha=.7)
                ax.set_axisbelow(True)
                ax.axvline(.5, color='.7', ls=':', lw=1)
                ax.set_ylabel(metric.capitalize() + ' (%)', fontsize=15)
                ax.set_xlabel('Distance to nearest reference PSV (bp; categorical bins)', fontsize=14, labelpad=12)
                scope_label = 'Segmental duplications' if scope == 'Segdup' else 'Segmental duplications ∩ low mappability'
                ax.set_title(f'{scope_label} | {vtype} | {metric.capitalize()}', fontsize=16, pad=42)
                ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.015), ncol=2,
                          frameon=False, fontsize=13, borderaxespad=0)
                fig.subplots_adjust(left=.095, right=.98, bottom=.19, top=.79)
                path = output / f'{scope}_{vtype}_{metric}_vs_psv_distance.png'
                fig.savefig(path, dpi=240, facecolor='white')
                files.append(str(path))
                plt.close(fig)
    grid = pd.DataFrame(plotted)
    # Every plotted point has the same category position across technologies.
    assert grid.groupby(['scope', 'vtype', 'metric', 'distance_bin']).x.nunique().eq(1).all()
    grid.to_csv(output.parent / 'plotted_values.tsv', sep='\t', index=False,
                na_rep='NaN', float_format='%.12g')
    (output.parent / 'plot_validation.json').write_text(json.dumps({
        'figures': files, 'n_figures': len(files), 'slots': len(grid),
        'finite_estimates': int(grid.value_percent.notna().sum()),
        'undefined_estimates': int(grid.value_percent.isna().sum()),
        'shared_categorical_x': True, 'confidence_intervals_drawn': False,
        'display_support_threshold': None, 'zero_denominator_is_missing': True,
        'zero_overlap_disconnected': True, 'y_ticks_percent': [0,20,40,60,80,100],
        'no_pdf_created': True}, indent=2) + '\n')
    return files
