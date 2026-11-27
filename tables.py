"""Tabular operations. Rows are features/taxa, columns are sample counts."""
import csv
import pandas as pd
import numpy as np

def load_counts(f):
    f.seek(0); first=f.readline(); f.seek(0)
    if isinstance(first,bytes): first=first.decode('utf-8-sig')
    sep='\t' if '\t' in first or first.startswith('# Constructed') else ','
    if first.startswith('# Constructed'):
        f.readline(); header=f.readline(); f.seek(0)
        if isinstance(header,bytes): header=header.decode('utf-8-sig')
    else: header=first
    cols=next(csv.reader([header],delimiter=sep))
    if len(cols)<2 or len(cols)!=len(set(cols)) or any(not c.strip() for c in cols[1:]): raise ValueError('Missing or duplicate sample columns.')
    d=pd.read_csv(f,sep=sep,skiprows=1 if first.startswith('# Constructed') else 0,dtype=str,encoding='utf-8-sig'); d=d.set_index(d.columns[0])
    if d.index.isna().any() or d.index.duplicated().any() or d.columns.duplicated().any(): raise ValueError('Feature and sample IDs must be unique and nonempty.')
    d.index=d.index.astype(str); d.columns=d.columns.astype(str); d=d.apply(pd.to_numeric,errors='raise')
    if not np.isfinite(d.to_numpy()).all() or (d<0).any().any(): raise ValueError('Counts must be finite and nonnegative.')
    return d

def merge_counts(tables,overlap='error'):
    if not tables: raise ValueError('Select tables')
    out=tables[0].copy()
    for d in tables[1:]:
        common=out.columns.intersection(d.columns)
        if len(common) and overlap=='error': raise ValueError('Overlapping sample IDs: '+', '.join(common[:10]))
        out=out.add(d,fill_value=0).fillna(0)
    return out

def collapse(d,tax,level):
    tax=tax.astype(str); missing=d.index.difference(tax.index)
    if len(missing): raise ValueError('Taxonomy missing '+str(len(missing))+' feature IDs.')
    def lineage(s):
        parts=[x.strip() for x in s.split(';')]; parts=[x if x and not x.endswith('__') else 'Unassigned' for x in parts]; return '; '.join((parts+['Unassigned']*7)[:level])
    return d.groupby(tax.reindex(d.index).map(lineage),sort=True).sum()

def aggregate(d,metadata,group,method):
    if metadata.index.duplicated().any(): raise ValueError('Duplicate metadata IDs')
    m=metadata.reindex(d.columns)[group]
    if m.isna().any() or (m.astype(str).str.strip()=='').any(): raise ValueError('Missing group assignments')
    return d.T.groupby(m).agg(method).T

def relative(d): return d.div(d.sum(axis=0).replace(0,np.nan),axis=1)*100
