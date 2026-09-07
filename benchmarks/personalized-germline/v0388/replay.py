#!/usr/bin/env python3
"""Replay discovery on exact selected lineage representatives, one animal at a time."""
import pathlib,subprocess,argparse,json
p=argparse.ArgumentParser();p.add_argument('--animal',choices=['ERR4238110','ERR4238104'],required=True);p.add_argument('--out',required=True);a=p.parse_args()
b=pathlib.Path(__file__).resolve().parent;repo=b.parents[2]
for arm in ['full','dropout','rare']:
 ref=b/('kimdb.V.fasta' if arm=='full' else f'{a.animal}.{arm}.V.fasta')
 out=pathlib.Path(a.out).resolve()/f'{a.animal}.{arm}'
 subprocess.run(['node',str(repo/'cli/swig-cli.mjs'),'personalized-germline','--airr',str(b/'inputs'/f'{a.animal}.{arm}.airr.tsv.gz'),'--v-reference',str(ref),'--out',str(out)],check=True)
 def active(path):
  d=json.loads(path.read_text());return sorted(x['sequence'] for p in d['pools'] for g in p['genes'] for x in g['activeAlleles'] if not x['known'])
 assert active(out/'personalized-germline.json')==active(b/'results'/f'{a.animal}.{arm}/personalized-germline.json'),f'Candidate mismatch: {arm}'
print('All three active novel sequence sets match the recorded results.')
