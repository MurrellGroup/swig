/** Global affine alignment: match 2, mismatch -3, gap opening -5, extension -1. */
export function globalCoordinateMap(left: string, right: string): { score: number; map: Int32Array } {
  const n=left.length,m=right.length,w=m+1,size=(n+1)*w;
  const M=new Float64Array(size).fill(-Infinity),X=M.slice(),Y=M.slice();
  const mt=new Uint8Array(size),xt=mt.slice(),yt=mt.slice();
  M[0]=0;
  for(let i=1;i<=n;i++) X[i*w]=-5-(i-1);
  for(let j=1;j<=m;j++) Y[j]=-5-(j-1);
  for(let i=1;i<=n;i++)for(let j=1;j<=m;j++){
    const k=i*w+j,d=k-w-1,u=k-w,l=k-1;
    const vals=[M[d],X[d],Y[d]];let best=0;for(let s=1;s<3;s++)if(vals[s]>vals[best])best=s;
    M[k]=vals[best]+(left[i-1]===right[j-1]?2:-3);mt[k]=best;
    if(M[u]-5>=X[u]-1){X[k]=M[u]-5;xt[k]=0;}else{X[k]=X[u]-1;xt[k]=1;}
    if(M[l]-5>=Y[l]-1){Y[k]=M[l]-5;yt[k]=0;}else{Y[k]=Y[l]-1;yt[k]=2;}
  }
  let i=n,j=m,k=n*w+m,state=0;const end=[M[k],X[k],Y[k]];
  for(let s=1;s<3;s++)if(end[s]>end[state])state=s;
  const score=end[state],map=new Int32Array(n).fill(-1);
  while(i||j){k=i*w+j;if(!i){j--;continue;}if(!j){i--;continue;}
    if(state===0){map[i-1]=j-1;state=mt[k];i--;j--;}else if(state===1){state=xt[k];i--;}else{state=yt[k];j--;}
  }
  return {score,map};
}
const eligibleCache=new Map<string,boolean>();
/** Ties permit an explicit ungapped optimal path; a strictly better indel path does not. */
export function substitutionCompatible(left:string,right:string,maximum=8):boolean{
  if(left.length!==right.length)return false;
  let differences=0;for(let i=0;i<left.length;i++)if(left[i]!==right[i]&&++differences>maximum)return false;
  if(!differences)return true;
  const key=left<right?left+'|'+right:right+'|'+left;
  const cached=eligibleCache.get(key);if(cached!==undefined)return cached;
  const eligible=2*left.length-5*differences>=globalCoordinateMap(left,right).score;
  eligibleCache.set(key,eligible);return eligible;
}
