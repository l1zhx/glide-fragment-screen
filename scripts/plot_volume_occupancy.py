#!/usr/bin/env python3
"""Annotate EXISTING ARM4-5 MDpocket volumes with MED occupancy.
No MDpocket execution, parameter changes, reference line, zero filling or p values.
Adapted from the project's 2026-10-27 arm45_volume_binding.py plotting flow.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

LIGANDS = ("AMP","BEN","CFM","CYC","CZX","DEX","DOP","ETI","HAB","HDQ","LEV","MED",
           "MEH","PRE","SAL","SOD","SUL","XYL")

def parse_descriptors(path):
    lines = [line.strip().lstrip("#").strip() for line in Path(path).read_text().splitlines()
             if line.strip()]
    header_index = next((i for i, line in enumerate(lines) if "pock_volume" in line.split()), None)
    if header_index is None:
        raise ValueError(f"Cannot find pock_volume header in {path}")
    header = lines[header_index].split()
    rows = []
    for line in lines[header_index + 1:]:
        tokens = line.split()
        if not tokens or tokens[0].startswith("#"):
            continue
        if len(tokens) != len(header):
            raise ValueError(f"Malformed descriptor row: {line}")
        rows.append(dict(zip(header, map(float, tokens))))
    if not rows:
        raise ValueError("Empty descriptor file")
    return rows


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path, rows, columns=None):
    columns = columns or (list(rows[0]) if rows else [])
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def new_folder(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def parse_bool(value):
    token=str(value).strip().lower()
    if token in {'1','1.0','true'}: return True
    if token in {'0','0.0','false'}: return False
    raise ValueError(f'Invalid occupancy value: {value!r}; cannot replace missing state by unbound')


def match_times(volume_times, state_times, tolerance=.001):
    state_times=np.asarray(state_times,float); volume_times=np.asarray(volume_times,float)
    if (not np.all(np.isfinite(state_times)) or not np.all(np.isfinite(volume_times))
            or len(state_times)<2 or np.any(np.diff(state_times)<=0)):
        raise ValueError('Occupancy times must be unique and strictly increasing')
    insertion=np.searchsorted(state_times,volume_times)
    left=np.clip(insertion-1,0,len(state_times)-1);right=np.clip(insertion,0,len(state_times)-1)
    idx=np.where(abs(state_times[left]-volume_times)<=abs(state_times[right]-volume_times),left,right)
    if np.any(abs(state_times[idx]-volume_times)>tolerance):
        raise ValueError('Occupancy and volume times do not match. Check replica, time units and time origin.')
    return idx


def state_spans(times, states):
    times=np.asarray(times,float); states=np.asarray(states,bool)
    if len(times)<2 or np.any(np.diff(times)<=0): raise ValueError('Invalid times')
    edges=np.r_[times[0]-(times[1]-times[0])/2,(times[:-1]+times[1:])/2,
                times[-1]+(times[-1]-times[-2])/2]
    changes=np.flatnonzero(np.diff(np.r_[False,states,False].astype(int)))
    return [(float(edges[a]),float(edges[b]-edges[a])) for a,b in zip(changes[::2],changes[1::2])]


def conditional_summary(v,detected,selected,label):
    v=np.asarray(v)[selected];d=np.asarray(detected)[selected]
    positive=v[d & (v>0)]
    return dict(group=label,n_frames=len(v),n_detected=int(d.sum()),
        detection_fraction=float(d.mean()) if len(d) else None,
        mean_raw_including_detection_zeros_A3=float(v.mean()) if len(v) else None,
        mean_detected_positive_A3=float(positive.mean()) if len(positive) else None,
        median_detected_positive_A3=float(np.median(positive)) if len(positive) else None,
        q25_detected_positive_A3=float(np.quantile(positive,.25)) if len(positive) else None,
        q75_detected_positive_A3=float(np.quantile(positive,.75)) if len(positive) else None)


def annotate(args):
    fields, volumes=read_csv(args.volume)
    occfields, occ=read_csv(args.occupancy)
    if not volumes or not occ: raise ValueError('Empty volume or occupancy file')
    key=args.state_column or f'{args.ligand}_ARM4-5'
    if key not in occfields: raise ValueError(f'Missing {key} in occupancy file')
    t=np.asarray([float(r['time_ns']) for r in volumes]); ot=np.asarray([float(r['time_ns']) for r in occ])
    if len(t)<2 or np.any(np.diff(t)<=0): raise ValueError('Volume times must increase')
    idx=match_times(t,ot,args.time_tolerance_ns)
    bound_full=np.array([parse_bool(r[key]) for r in occ]);bound=bound_full[idx]
    other_keys=[f'{name}_ARM4-5' for name in LIGANDS if name!=args.ligand]
    all_others_known=all(k in occfields for k in other_keys)
    other_full=np.any([[parse_bool(r[k]) for r in occ] for k in other_keys],axis=0) if all_others_known else None
    other=other_full[idx] if other_full is not None else None
    v=np.array([float(r['volume_A3']) for r in volumes])
    if 'nb_AS' in fields:
        sphere_counts=np.array([float(r['nb_AS']) for r in volumes])
        if not np.all(np.isfinite(sphere_counts)) or np.any(sphere_counts<0):
            raise ValueError('Invalid alpha-sphere counts')
        detected=sphere_counts>0;detection_source='nb_AS'
    else:
        descriptor=getattr(args,'descriptors',None) or args.volume.parent/'measure/mdpout_descriptors.txt'
        if not descriptor.is_file(): raise FileNotFoundError('Need nb_AS column or adjacent measure/mdpout_descriptors.txt')
        records=parse_descriptors(descriptor)
        if len(records)!=len(v)+1 or not np.array_equal([r['pock_volume'] for r in records[1:]],v):
            raise ValueError('Descriptor/volume mismatch')
        detected=np.array([r['nb_AS']>0 for r in records[1:]]);detection_source=str(descriptor)
    if not np.all(np.isfinite(v)) or np.any(v<0):raise ValueError('Invalid volume')
    if np.any((~detected) & (v>0)):
        raise ValueError('Positive volume with zero alpha spheres: verify input correspondence')
    out=new_folder(args.out)
    joined=[]
    for i,row in enumerate(volumes):
        joined.append(dict(**row,MED_contact_defined=int(bound[i]),
            other_fragments_contact_defined=int(other[i]) if other is not None else 'unknown',
            cavity_detected=int(detected[i]),occupancy_time_ns=ot[idx[i]]))
    write_csv(out/'volume_occupancy_joined.csv',joined)
    groups=[('MED present',bound),('MED absent',~bound)]
    if other is not None:
        groups += [('MED only',bound & ~other),('MED and others',bound & other),
                   ('Other fragments only',~bound & other),('No fragment detected',~bound & ~other)]
    stats=[conditional_summary(v,detected,s,label) for label,s in groups]
    write_csv(out/'conditional_volume_summary.csv',stats)
    def difference(key):
        x,y=stats[0][key],stats[1][key]
        return x-y if x is not None and y is not None else None
    save_json(out/'MED_present_minus_absent.json',dict(
        delta_detection_fraction=difference('detection_fraction'),
        delta_detected_positive_mean_A3=difference('mean_detected_positive_A3'),
        delta_detected_positive_median_A3=difference('median_detected_positive_A3'),
        interpretation='Descriptive association only; detection selection, time drift, other fragments and autocorrelation may confound.'))
    intervals=[dict(start_ns=a,end_ns=a+b,state='MED contact-defined') for a,b in state_spans(ot,bound_full)]
    write_csv(out/'MED_intervals.csv',intervals,['start_ns','end_ns','state'])
    save_json(out/'analysis_provenance.json',dict(volume=str(args.volume.resolve()),volume_sha256=sha(args.volume),
        occupancy=str(args.occupancy.resolve()),occupancy_sha256=sha(args.occupancy),state_column=key,
        time_tolerance_ns=args.time_tolerance_ns,detection_source=detection_source,
        other_fragment_states_complete=all_others_known,ligand=args.ligand,
        caveats=['Contact-defined occupancy is not affinity.',
                 'MED absent in mixed-fragment simulation is not an apo control.',
                 'Frames are time-correlated; no independent-frame p values.',
                 'Early binding and late absence confound binding state with time.',
                 'All-zero volumes preserved; plotted as no matched cavity, not forced closure.']))
    plots(out,t,v,detected,ot,bound_full,other_full,bound,stats,args.ligand)
    print(f'Annotated {len(v)} sampled frames; MED present in {bound.sum()}. Others known: {all_others_known}',flush=True)


def plots(out,t,v,detected,ot,bound_full,other_full,bound,stats,ligand,reference_volume=None):
    # reference_volume is accepted for compatibility, but not displayed.
    # Source reference records and measurement metadata remain unchanged.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,(ax,stateax)=plt.subplots(2,1,figsize=(11,5),sharex=True,
        gridspec_kw={'height_ratios':[4,1]},constrained_layout=True)
    for start,width in state_spans(ot,bound_full):
        ax.axvspan(start/1000,(start+width)/1000,color='#D55E00',alpha=.12,lw=0)
    good=detected & (v>0)
    ax.plot(t/1000,np.where(good,v,np.nan),color='#167D8D',lw=.8,marker='.',markersize=2,
            label='Detected cavity volume')
    if np.any(~good):ax.scatter(t[~good]/1000,np.zeros((~good).sum()),s=14,marker='x',color='#777777',label='No matched cavity / zero estimate')
    ax.set(ylabel=r'MDpocket volume ($\mathrm{\AA}^3$)',title='ARM4-5: volume and contact-defined occupancy')
    ax.legend(frameon=False,loc='upper right')
    ax.text(.01,.97,f'{ligand} occupancy: shaded intervals; zero estimates retained',
            transform=ax.transAxes,va='top',fontsize=9)
    bands=[(ligand,bound_full,'#D55E00')]
    if other_full is not None:bands.append(('Any other fragment',other_full,'#7B6BA8'))
    for row,(label,state,color) in enumerate(bands):
        spans=[(a/1000,b/1000) for a,b in state_spans(ot,state)]
        stateax.broken_barh(spans,(row-.3,.6),facecolors=color)
    stateax.set(yticks=list(range(len(bands))),yticklabels=[x[0] for x in bands],
                xlabel=r'Time ($\mu$s)',xlim=(t[0]/1000,t[-1]/1000),ylim=(-.6,len(bands)-.4))
    for ext in ('png','svg'):fig.savefig(out/f'volume_with_MED_occupancy.{ext}',dpi=350)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
    labels=[f'{ligand} present',f'{ligand} absent']
    for i,r in enumerate(stats[:2]):
        if r['detection_fraction'] is not None:
            axes[0].bar(i,r['detection_fraction'],color=['#D55E00','#167D8D'][i])
        else:
            axes[0].text(i,.5,'N/A: no frames',ha='center',rotation=90)
    axes[0].set(xticks=[0,1],xticklabels=labels,ylim=(0,1.15),ylabel='Cavity detection fraction')
    for i,r in enumerate(stats[:2]):axes[0].text(i,1.03,f'{r["n_detected"]}/{r["n_frames"]}',ha='center')
    for pos,sel in enumerate([bound,~bound],1):
        data=v[sel & good]
        if len(data):axes[1].boxplot([data],positions=[pos],widths=.5,showfliers=False)
        else:axes[1].text(pos,0,'No detected frames',ha='center',rotation=90)
    axes[1].set(xticks=[1,2],xticklabels=labels,ylabel=r'Detected-positive volume ($\mathrm{\AA}^3$)',
                title='Conditional distribution; correlated frames')
    for ext in ('png','svg'):fig.savefig(out/f'MED_present_absent_comparison.{ext}',dpi=350)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--volume',type=Path,required=True)
    p.add_argument('--occupancy',type=Path,required=True)
    p.add_argument('--descriptors',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--ligand',choices=['MED'],default='MED')
    p.add_argument('--state-column',default='MED_ARM4-5')
    p.add_argument('--time-tolerance-ns',type=float,default=.001)
    a=p.parse_args()
    if not np.isfinite(a.time_tolerance_ns) or a.time_tolerance_ns<=0:
        p.error('Time tolerance must be finite and positive')
    annotate(a)

if __name__=='__main__':
    main()

