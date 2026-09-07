#!/usr/bin/env python3
import json,pathlib,math,sys
b=pathlib.Path(__file__).resolve().parent;old=b.parent/'v0388';animal=sys.argv[1]
def read(arm):return json.loads((b/f'{animal}.{arm}/joint-germline.json').read_text())
def evidence(d):
 out={}
 for g in d['results']:
  threshold=math.log(max(1,d['evidenceRule']['testsBySubject'][g['key'].split('|')[0]])/d['evidenceRule']['alpha'])
  for e in g['evidence']:
   out[e['sequence']]={**e,'logThreshold':threshold,'positive':e['logEvidence'] is not None and e['logEvidence']>0,'passes':e['logEvidence'] is not None and e['logEvidence']>=threshold}
 return out
full=read('full');fe=evidence(full);result={'animal':animal,'full':{'positive':sum(e['positive'] for e in fe.values()),'passes':sum(e['passes'] for e in fe.values())},'panels':{}}
for arm in ['dropout','rare']:
 d=read(arm);es=evidence(d);panel=json.loads((old/f'{animal}{".rare" if arm=="rare" else ""}.panel.json').read_text())['panel'];targets={x['sequence'][:-2] for x in panel};rows=[]
 for x in panel:
  e=es.get(x['sequence'][:-2]);rows.append({'target':x['allele'],'nearExactCdr3s':x['near_cdr3'],'nominated':e is not None,'logEvidence':e['logEvidence'] if e else None,'positive':e['positive'] if e else False,'passes':e['passes'] if e else False,'nullOptimizationGap':e.get('gap') if e else None,'logThreshold':e.get('logThreshold') if e else None})
 extras=[e for seq,e in es.items() if seq not in targets and e['positive']]
 result['panels'][arm]={'targets':rows,'nominated':sum(r['nominated'] for r in rows),'positive':sum(r['positive'] for r in rows),'passes':sum(r['passes'] for r in rows),'positiveExtras':extras,'passingExtras':[e for e in extras if e['passes']],'positiveExtrasAbsentFull':[e for e in extras if e['sequence'] not in fe or not fe[e['sequence']]['positive']]}
(b/f'{animal}.score.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
