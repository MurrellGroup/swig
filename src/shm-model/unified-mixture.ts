/** One normalized model, one objective. Experimental; not a production caller. */
export interface Row { a:Float64Array; b:Float64Array; count?:number }
export interface Fit { theta:Float64Array;q:Float64Array;ll:number;upper:number;gap:number;iterations:number }
export function likelihood(rows:readonly Row[],theta:ArrayLike<number>,q:ArrayLike<number>):number {
 let ll=0;for(const r of rows){let p=0;for(let c=0;c<theta.length;c++)p+=theta[c]*r.a[c]+q[c]*r.b[c];ll+=(r.count??1)*Math.log(Math.max(1e-300,p));}return ll;
}
/** Product-of-simplexes concave MLE; tangent-plane bound certifies an upper bound
 * even if optimization stops early. No local optimum can exaggerate evidence. */
export function fit(rows:readonly Row[],allowed:readonly number[],initial?:Fit,maxIterations=1000,tolerance=1e-5):Fit {
 const C=rows[0].a.length,theta=new Float64Array(C),q=new Float64Array(C),n=rows.reduce((s,r)=>s+(r.count??1),0);
 for(const c of allowed)theta[c]=Math.max(1e-6,initial?.theta[c]??1/allowed.length);
 for(let c=0;c<C;c++)q[c]=Math.max(1e-6,initial?.q[c]??1/C);
 const normalize=(x:Float64Array)=>{const s=x.reduce((a,b)=>a+b,0);for(let i=0;i<x.length;i++)x[i]/=s;};normalize(theta);normalize(q);
 let ll=-Infinity,gap=Infinity,iteration=0;
 for(;iteration<maxIterations;iteration++){
  const ga=new Float64Array(C),gb=new Float64Array(C);ll=0;
  for(const r of rows){let p=0;for(let c=0;c<C;c++)p+=theta[c]*r.a[c]+q[c]*r.b[c];p=Math.max(1e-300,p);const w=(r.count??1)/p;ll+=(r.count??1)*Math.log(p);for(let c=0;c<C;c++){ga[c]+=w*r.a[c];gb[c]+=w*r.b[c];}}
  let da=0,db=0,ma=-Infinity,mb=-Infinity;for(const c of allowed){da+=theta[c]*ga[c];ma=Math.max(ma,ga[c]);}for(let c=0;c<C;c++){db+=q[c]*gb[c];mb=Math.max(mb,gb[c]);}
  gap=Math.max(0,ma-da)+Math.max(0,mb-db);
  if(gap<tolerance||iteration===maxIterations-1)break;
  // EM preserves positivity; revive tiny coordinates to avoid a false boundary
  // convergence. The bound above includes every allowed coordinate regardless.
  for(const c of allowed)theta[c]=Math.max(1e-15,theta[c]*ga[c]/Math.max(1e-300,da));
  for(let c=0;c<C;c++)q[c]=Math.max(1e-15,q[c]*gb[c]/Math.max(1e-300,db));normalize(theta);normalize(q);
 }
 return {theta,q,ll,upper:ll+gap,gap,iterations:iteration+1};
}
export const HAZARDS=[0,1,10,100,Infinity] as const;
/** logKernel is flattened candidate × exposure. A and B use a single row scale
 * shared across ALL hazard models, essential for valid likelihood comparisons. */
export function integrate(logKernel:Float64Array,times:readonly number[],weights:readonly number[]):Row[]{
 const K=times.length,C=logKernel.length/K,max=Math.max(...logKernel),kernel=Float64Array.from(logKernel,v=>Math.exp(v-max));
 return HAZARDS.map(h=>{
  const a=new Float64Array(C),b=new Float64Array(C);
  for(let c=0;c<C;c++)for(let k=0;k<K;k++){
   const retained=times[k]===0?1:h===Infinity?0:Math.exp(-h*times[k]);const value=weights[k]*kernel[c*K+k];a[c]+=retained*value;b[c]+=(1-retained)*value;
  }return {a,b};
 });
}
export function fitHazards(rows:readonly Row[][],allowed:readonly number[],initial?:Fit,maxIterations=1000){
 const fits=rows.map(r=>fit(r,allowed,initial,maxIterations));let best=0;for(let h=1;h<fits.length;h++)if(fits[h].ll>fits[best].ll)best=h;
 return {best,fit:fits[best],upper:Math.max(...fits.map(f=>f.upper)),fits};
}
export function splitEvidence(train:readonly Row[][],test:readonly Row[][],known:readonly boolean[],maxIterations=1000){
 const all=known.map((_,i)=>i),alternative=fitHazards(train,all,undefined,maxIterations),predictive=likelihood(test[alternative.best],alternative.fit.theta,alternative.fit.q);
 const evidence=[];
 for(let c=0;c<known.length;c++)if(!known[c]){
  // No test-data-dependent nomination: candidates with zero fitted inheritance
  // need no expensive test and receive an evidence value of zero.
  if(alternative.fit.theta[c]<1e-6){evidence.push({candidate:c,logEvidence:-Infinity,nullUpper:null,gap:null});continue;}
  // Average two pre-fitted predictors: a broad allele model and a focused
  // known-plus-candidate model. A mixture of normalized predictors is normalized;
  // unlike choosing the better test-set score, it pays its mixture weight.
  const focused=fit(train[alternative.best],all.filter(i=>known[i]||i===c),alternative.fit,maxIterations);
  const focusedPredictive=likelihood(test[alternative.best],focused.theta,focused.q);
  const m=Math.max(predictive,focusedPredictive),candidatePredictive=m+Math.log((Math.exp(predictive-m)+Math.exp(focusedPredictive-m))/2);
  const allowed=all.filter(i=>i!==c),feasibleTheta=Float64Array.from(alternative.fit.theta);feasibleTheta[c]=0;
  const sum=feasibleTheta.reduce((a,b)=>a+b,0);for(const i of allowed)feasibleTheta[i]=sum>0?feasibleTheta[i]/sum:1/allowed.length;
  const feasible=Math.max(...test.map(rows=>likelihood(rows,feasibleTheta,alternative.fit.q)));
  // A feasible null already beating the predictor proves positive evidence is
  // impossible. Report zero evidence rather than an unbounded optimization score.
  if(feasible>=candidatePredictive){evidence.push({candidate:c,logEvidence:-Infinity,nullUpper:null,gap:null});continue;}
  const nullFit=fitHazards(test,allowed,alternative.fit,maxIterations);
  evidence.push({candidate:c,logEvidence:candidatePredictive-nullFit.upper,nullUpper:nullFit.upper,gap:Math.max(...nullFit.fits.map(f=>f.gap))});
 }
 return {hazard:HAZARDS[alternative.best],trainingLogLikelihood:alternative.fit.ll,trainingGap:alternative.fit.gap,predictiveLogLikelihood:predictive,theta:Array.from(alternative.fit.theta),somatic:Array.from(alternative.fit.q),evidence};
}
