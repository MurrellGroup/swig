import { createReadStream } from 'node:fs';
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { createGunzip } from 'node:zlib';
import { createInterface } from 'node:readline';
import { join } from 'node:path';
import { PersonalizedGermlineAccumulator, personalizedGermlineEvidenceRows, personalizedGermlineFasta } from '../src/personalized-germline.ts';
export async function runPersonalizedGermline(args){
 const options={};for(let i=0;i<args.length;i++){if(args[i]==='--help'){console.log('swig-cli personalized-germline --airr FILE[.gz] --v-reference V.fasta --out DIRECTORY [--subject ID]\nRequires clone_id (or lineage_id), V alignments and germline coordinates. Raw sequence and v_sequence_end support boundary uncertainty. One model per subject.');return;}if(!['--airr','--v-reference','--out','--subject'].includes(args[i])||!args[i+1])throw new Error(`Unknown or incomplete personalized-germline option: ${args[i]}`);options[args[i]]=args[++i];}
 for(const key of ['--airr','--v-reference','--out'])if(!options[key])throw new Error(`personalized-germline requires ${key}.`);
 const fasta=await readFile(options['--v-reference'],'utf8'),acc=new PersonalizedGermlineAccumulator(fasta);
 const input=createReadStream(options['--airr']);const stream=options['--airr'].endsWith('.gz')?input.pipe(createGunzip()):input;
 input.on('error',error=>stream.destroy(error));
 let headers,ordinal=0;
 for await(const line of createInterface({input:stream,crlfDelay:Infinity})){
  if(!headers){headers=line.replace(/^\uFEFF/,'').split('\t');for(const name of ['v_call','v_sequence_alignment','v_germline_alignment','v_germline_start'])if(!headers.includes(name))throw new Error(`Missing AIRR field: ${name}`);if(!headers.includes('clone_id')&&!headers.includes('lineage_id'))throw new Error('Input must contain clone_id or lineage_id from lineage assignment.');continue;}
  if(!line)continue;const cells=line.split('\t'),row=Object.fromEntries(headers.map((h,i)=>[h,cells[i]??'']));if(options['--subject'])row.subject_id=options['--subject'];acc.add(row,ordinal++,Number(row.clone_id||row.lineage_id));
 }
 const start=performance.now();let last=0;
 const result=acc.finish((done,total)=>{if(performance.now()-last>5000){process.stderr.write(`[personalized] ${done}/${total} competitive components\n`);last=performance.now();}});
 const seconds=(performance.now()-start)/1000;await mkdir(options['--out'],{recursive:true});
 await writeFile(join(options['--out'],'personalized-germline.json'),JSON.stringify({...result,fitSeconds:seconds},null,2));
 const evidence=personalizedGermlineEvidenceRows(result),keys=Object.keys(evidence[0]??{});
 await writeFile(join(options['--out'],'personalized-germline.tsv'),keys.join('\t')+'\n'+evidence.map(row=>keys.map(k=>String(row[k]??'').replace(/[\t\r\n]/g,' ')).join('\t')).join('\n')+'\n');
 for(const pool of result.pools)await writeFile(join(options['--out'],`${pool.id}.V.fasta`),personalizedGermlineFasta(fasta,result,pool.id));
 console.log(JSON.stringify({representativeLineages:result.representativeLineages,novel:result.activeNovelAlleles,known:result.activeKnownAlleles,fitSeconds:seconds,warnings:result.warnings},null,2));
}
