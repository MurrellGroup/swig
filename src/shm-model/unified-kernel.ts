import {hs5fRate} from './hs5f.ts';
const I=()=>Float64Array.from({length:16},(_,i)=>Number(i%5===0));
/** Exact (to Poisson truncation tolerance) four-base CTMC, frozen flanking context.
 * Alternative bases are states of the same process; returns and repeated hits
 * are possible. Sequencing error is a separate stochastic channel. */
export function transition(context:string,t:number,error=.001):Float64Array {
 const Q=new Float64Array(16);let rate=0;
 for(let a=0;a<4;a++){const seq=context.slice(0,2)+'ACGT'[a]+context.slice(3);let sum=0;for(let b=0;b<4;b++)if(a!==b){Q[a*4+b]=hs5fRate(seq,2,'ACGT'[b]);sum+=Q[a*4+b];}Q[a*4+a]=-sum;rate=Math.max(rate,sum);}
 const P=I();if(t>0&&rate>0){
  const R=Float64Array.from(Q,(v,i)=>v/rate+Number(i%5===0));let power=I(),poisson=Math.exp(-rate*t),mass=poisson;for(let i=0;i<16;i++)P[i]=power[i]*poisson;
  for(let n=1;n<1000&&1-mass>1e-14;n++){
   const next=new Float64Array(16);for(let a=0;a<4;a++)for(let b=0;b<4;b++)for(let k=0;k<4;k++)next[a*4+b]+=power[a*4+k]*R[k*4+b];power=next;poisson*=rate*t/n;mass+=poisson;for(let i=0;i<16;i++)P[i]+=poisson*power[i];
  }
  if(mass<1-1e-10)throw new Error('CTMC quadrature did not converge');
  for(let a=0;a<4;a++){let s=0;for(let b=0;b<4;b++)s+=P[a*4+b];for(let b=0;b<4;b++)P[a*4+b]/=s;}
 }
 return Float64Array.from(P,p=>error/3+(1-4*error/3)*p);
}
