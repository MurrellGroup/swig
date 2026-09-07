#!/usr/bin/env python3
"""Rebuild the exact three arms from a reuploaded FASTA, using fixed recorded panels."""
import pathlib,argparse,gzip,random,json,subprocess,hashlib
p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--animal',required=True,choices=['ERR4238110','ERR4238104']);p.add_argument('--out',required=True);args=p.parse_args()
b=pathlib.Path(__file__).resolve().parent;repo=b.parents[2];out=pathlib.Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
rng=random.Random(1);reservoir=[];n=0
source=pathlib.Path(args.fasta)
def records(f):
 header=None;seq=[]
 for line in f:
  if line.startswith('>'):
   if header is not None:yield header+''.join(seq)
   header=line;seq=[]
  else:seq.append(line)
 if header is not None:yield header+''.join(seq)
with (gzip.open(source,'rt') if source.suffix=='.gz' else source.open()) as f:
 for n,record in enumerate(records(f),1):
  item=(n-1,record)
  if n<=300000:reservoir.append(item)
  else:
   j=rng.randrange(n)
   if j<300000:reservoir[j]=item
sample=out/f'{args.animal}.300k.fasta';sample.write_text(''.join(record for _,record in sorted(reservoir)))
expected={'ERR4238110':'6a87b2a2caf75c5bfbeb5620210afab0b2a9a781d0c07a32b9e0c370f297af8a','ERR4238104':'287e74867e84f5238a8f9d348a454ff772ae815b038f88e16fa9eba1780ada62'}
assert hashlib.sha256(sample.read_bytes()).hexdigest()==expected[args.animal], 'Input/sample differs from recorded validation'
for arm in ['full','dropout','rare']:
 cfg=json.loads((b/f'{args.animal}.{arm}.config.json').read_text());cfg['inputs'][0]['path']=str(sample)
 cfg['references']['files']={s:str(b/name) for s,name in cfg['references']['files'].items()};cfg['output']['directory']=str(out/f'{args.animal}.{arm}')
 config=out/f'{args.animal}.{arm}.config.json';config.write_text(json.dumps(cfg,indent=2))
 cli=['node',str(repo/'cli/swig-cli.mjs')]
 subprocess.run(cli+['run','--config',str(config)],check=True)
 subprocess.run(cli+['personalized-germline','--airr',str(pathlib.Path(cfg['output']['directory'])/f'{args.animal}.processed.airr.tsv.gz'),'--v-reference',cfg['references']['files']['V'],'--out',str(pathlib.Path(cfg['output']['directory'])/'discovery')],check=True)
