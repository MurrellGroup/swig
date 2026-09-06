#!/usr/bin/env python3
"""Score exact recovery, terminal-only differences, and paired-control extras."""
import json,pathlib,sys
r=pathlib.Path(sys.argv[1]);animal=sys.argv[2]
def load(mode):return json.loads((r/f'{animal}.{mode}'/'validated/personalized-germline.json').read_text())
full,drop=load('full'),load('dropout')
def novel(d):return {a['sequence']:a for pool in d['pools'] for g in pool['genes'] for a in g['activeAlleles'] if not a['known']}
f,d=novel(full),novel(drop);panel=json.loads((r/f'{animal}.panel.json').read_text())['panel'];recoveries=[];selected=set()
for target in panel:
 candidates=sorted([(sum(a!=b for a,b in zip(sequence,target['sequence'])),sequence,allele) for sequence,allele in d.items() if len(sequence)==len(target['sequence'])],key=lambda item:(item[0],item[2]['id']))
 if not candidates:recoveries.append({'target':target['allele'],'status':'miss'});continue
 distance,sequence,allele=candidates[0];differences=[i+1 for i,(a,b) in enumerate(zip(sequence,target['sequence'])) if a!=b]
 status='exact' if distance==0 else 'terminal-only' if all(p>len(sequence)-12 for p in differences) else 'miss'
 if status!='miss':selected.add(sequence)
 recoveries.append({'target':target['allele'],'status':status,'nearestCandidate':allele['id'],'differences':differences,'nearExactBaselineCdr3s':target['near_cdr3'],'nearestSurvivingDistance':target['distance']})
extras=[allele['id'] for sequence,allele in d.items() if sequence not in f and sequence not in selected]
result={'animal':animal,'recoveries':recoveries,'exact':sum(x['status']=='exact' for x in recoveries),'terminalOnly':sum(x['status']=='terminal-only' for x in recoveries),'dropoutInducedExtras':extras,'sharedNonTargetCandidates':[allele['id'] for sequence,allele in d.items() if sequence in f and sequence not in selected],'fullOnlyCandidates':[allele['id'] for sequence,allele in f.items() if sequence not in d],'fullNovel':len(f),'dropoutNovel':len(d),'fullLineages':full['representativeLineages'],'dropoutLineages':drop['representativeLineages'],'fullFitSeconds':full['fitSeconds'],'dropoutFitSeconds':drop['fitSeconds']}
(r/f'{animal}.score.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
