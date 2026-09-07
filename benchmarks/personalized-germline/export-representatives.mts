import {createReadStream,createWriteStream,readFileSync} from 'node:fs';
import {createGunzip,createGzip} from 'node:zlib';
import {createInterface} from 'node:readline';
import {once} from 'node:events';
import {pipeline} from 'node:stream/promises';
import {PersonalizedGermlineAccumulator} from '../../src/personalized-germline.ts';
const [source,reference,destination]=process.argv.slice(2);
if(!source||!reference||!destination)throw new Error('Usage: export-representatives.mts processed.airr.tsv.gz V.fasta representatives.airr.tsv.gz');
const fields=['sequence_id','subject_id','clone_id','lineage_id','v_call','v_sequence_alignment','v_germline_alignment','v_germline_start','v_sequence_start','v_sequence_end','sequence','rev_comp','d_sequence_start','j_sequence_start'];
const acc=new PersonalizedGermlineAccumulator(readFileSync(reference,'utf8'));
async function* rows(){
 let headers:string[]=[];
 const stream=createReadStream(source).pipe(createGunzip());
 for await(const line of createInterface({input:stream,crlfDelay:Infinity})){
  if(!headers.length){headers=line.split('\t');continue;}if(!line)continue;
  const cells=line.split('\t');yield Object.fromEntries(fields.map(f=>[f,cells[headers.lastIndexOf(f)]??'']));
 }
}
let ordinal=0;for await(const row of rows())acc.add(row,ordinal++,Number(row.clone_id||row.lineage_id));
const selected=new Set(acc.selectedRepresentativeOrdinals()),gzip=createGzip(),done=pipeline(gzip,createWriteStream(destination));
gzip.write(fields.join('\t')+'\n');ordinal=0;
for await(const row of rows())if(selected.has(ordinal++))if(!gzip.write(fields.map(f=>row[f]).join('\t')+'\n'))await once(gzip,'drain');
gzip.end();await done;console.log(JSON.stringify({source,representatives:selected.size,destination}));
