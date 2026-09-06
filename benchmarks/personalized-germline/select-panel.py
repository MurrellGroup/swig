import csv,gzip,json,sys,pathlib,collections
r=pathlib.Path(__file__).parent;a=sys.argv[1]
refs={};head=None
for line in (r/'kimdb.V.fasta').read_text().splitlines():
 if line.startswith('>'):head=line[1:];refs[head]=[]
 else:refs[head].append(line)
refs={h:''.join(s) for h,s in refs.items()};byname={h.split()[0]:(h,s) for h,s in refs.items()};copies=collections.Counter(refs.values());stats=collections.defaultdict(lambda:{'reads':0,'cdr3':set(),'near':set()})
with gzip.open(r/f'{a}.full'/f'{a}.annotated.airr.tsv.gz','rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  v=row['v_call'];cdr3=row.get('cdr3') or row.get('junction');q=row['v_sequence_alignment'];g=row['v_germline_alignment'];pairs=[(x,y) for x,y in zip(q,g) if x in 'ACGT' and y in 'ACGT'];shm=sum(x!=y for x,y in pairs)/max(1,len(pairs))
  if v not in byname or not cdr3 or len(pairs)<240:continue
  st=stats[v];st['reads']+=1;st['cdr3'].add(cdr3)
  if shm<=.01:st['near'].add(cdr3)
candidates=[]
for name,st in stats.items():
 seq=byname[name][1]
 if copies[seq]!=1 or len(st['near'])<20:continue
 distances=[(sum(x!=y for x,y in zip(seq,s)),other) for other,(h,s) in byname.items() if other!=name and other.split('*')[0]==name.split('*')[0] and len(seq)==len(s)]
 if not distances:continue
 distance,neighbor=min(distances)
 if 1<=distance<=6:candidates.append({'allele':name,'nearest':neighbor,'distance':distance,'reads':st['reads'],'near_cdr3':len(st['near']),'cdr3':len(st['cdr3']),'sequence':seq})
candidates.sort(key=lambda x:(-x['near_cdr3'],x['allele']))
panel=[];genes=set()
for distance in [6,5,4,3,2,1]:
 eligible=[x for x in candidates if x['distance']==distance and x['allele'].split('*')[0] not in genes]
 if eligible:panel.append(eligible[0]);genes.add(eligible[0]['allele'].split('*')[0])
for x in candidates:
 if len(panel)>=6:break
 if x['allele'].split('*')[0] not in genes:panel.append(x);genes.add(x['allele'].split('*')[0])
assert len(panel)==6
(r/f'{a}.panel.json').write_text(json.dumps({'selection':'Sequence-unique; >=20 distinct near-exact CDR3s (<=1% V SHM); >=240 aligned bases; same-gene same-length surviving neighbor within 1-6 SNPs. One per gene, strongest support within distance strata; fill remaining by support.','candidates':candidates,'panel':panel},indent=2))
dropped={x['allele'] for x in panel};(r/f'{a}.dropout.V.fasta').write_text(''.join(f'>{h}\n{s}\n' for h,s in refs.items() if h.split()[0] not in dropped))
config=json.loads((r/f'{a}.full.config.json').read_text());config['references']['files']['V']=str((r/f'{a}.dropout.V.fasta').resolve());config['output']['directory']=str((r/f'{a}.dropout').resolve());(r/f'{a}.dropout.config.json').write_text(json.dumps(config,indent=2))
print(json.dumps([{k:v for k,v in x.items() if k!='sequence'} for x in panel],indent=2))
