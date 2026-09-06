#!/usr/bin/env python3
"""Run one animal at a time; create paired 300k-read KIMDB dropout controls."""
import argparse,gzip,hashlib,json,pathlib,random,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--animal',required=True);p.add_argument('--workdir',required=True);p.add_argument('--workers',type=int,default=8);args=p.parse_args()
repo=pathlib.Path(__file__).resolve().parents[2];work=pathlib.Path(args.workdir).resolve();work.mkdir(parents=True,exist_ok=True);animal=args.animal
cli=['node',str(repo/'cli/swig-cli.mjs')]
def run(cmd):subprocess.run(cmd,check=True)
source=pathlib.Path(args.fasta).resolve();rng=random.Random(1);sample=[];count=0
with (gzip.open(source,'rt') if source.suffix=='.gz' else source.open()) as f:
 header=None;lines=[]
 def add(header,lines):
  global count
  item=(count,header+''.join(lines))
  if count<300000:sample.append(item)
  else:
   j=rng.randrange(count+1)
   if j<300000:sample[j]=item
  count+=1
 for line in f:
  if line.startswith('>'):
   if header is not None:add(header,lines)
   header=line;lines=[]
  else:lines.append(line)
 if header is not None:add(header,lines)
sample.sort();samplepath=work/f'{animal}.300k.fasta'
with samplepath.open('w') as f:
 for _,record in sample:f.write(record)
(work/f'{animal}.sampling.json').write_text(json.dumps({'input':str(source),'inputRecords':count,'selected':len(sample),'algorithm':'Python random.Random(1), Algorithm R, restore original ordinal order','sha256':hashlib.sha256(samplepath.read_bytes()).hexdigest()},indent=2))
ref=repo/'public/references/kimdb-1.1/Macaca_mulatta/IGH'
run(cli+['prepare-reference','-germline_db_V',str(ref/'V.fasta'),'-germline_db_D',str(ref/'D.fasta'),'-germline_db_J',str(ref/'J.fasta'),'-organism','rhesus_monkey','--out-prefix',str(work/'kimdb')])
config={'inputs':[{'path':str(samplepath),'sampleId':animal,'subjectId':animal}], 'references':{'species':'Macaca mulatta_AG07107','scope':'IGH','prepareMetadata':False,'files':{s:str(work/f'kimdb.{s}.fasta') for s in 'VDJ'}},'annotation':{'assignerStrategy':'aer_robust','callingProfile':'r_optimized','workers':args.workers},'preprocessing':{'subsample':{'enabled':False},'fastqFilter':{'enabled':False}},'pipeline':{'collapse':{'enabled':True,'mode':'indel'},'lineage':{'enabled':True},'shm':{'enabled':False}},'output':{'directory':str(work/f'{animal}.full'),'airrCompression':'gzip','prefix':animal}}
full=work/f'{animal}.full.config.json';full.write_text(json.dumps(config,indent=2));run(cli+['run','--config',str(full)])
# Selection reads only the full-reference assignment evidence, before discovery.
selector=pathlib.Path(__file__).with_name('select-panel.py').read_text();local=work/'select-panel.py';local.write_text(selector);run([sys.executable,str(local),animal])
run(cli+['run','--config',str(work/f'{animal}.dropout.config.json')])
for mode in ['full','dropout']:
 run(cli+['personalized-germline','--airr',str(work/f'{animal}.{mode}'/f'{animal}.processed.airr.tsv.gz'),'--v-reference',str(work/('kimdb.V.fasta' if mode=='full' else f'{animal}.dropout.V.fasta')),'--out',str(work/f'{animal}.{mode}'/'validated')])
print('Completed paired analysis for',animal)
