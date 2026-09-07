#!/usr/bin/env python3
import pathlib,subprocess,argparse
p=argparse.ArgumentParser();p.add_argument('--animal',required=True,choices=['ERR4238110','ERR4238104']);p.add_argument('--out',required=True);args=p.parse_args();b=pathlib.Path(__file__).resolve().parent;repo=b.parents[2];old=b.parent/'v0388'
for arm in ['full','dropout','rare']:
 ref=old/('kimdb.V.fasta' if arm=='full' else f'{args.animal}.{arm}.V.fasta')
 subprocess.run(['node',str(repo/'cli/swig-cli.mjs'),'joint-germline','--airr',str(old/'inputs'/f'{args.animal}.{arm}.airr.tsv.gz'),'--v-reference',str(ref),'--out',str(pathlib.Path(args.out).resolve()/f'{args.animal}.{arm}')],check=True)
