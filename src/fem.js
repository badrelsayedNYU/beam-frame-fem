/* FrameLab: deterministic linear elastic Euler–Bernoulli frame FEM. */
(function(root){
'use strict';
const zeros=(r,c)=>Array.from({length:r},()=>Array(c).fill(0));
const transpose=a=>a[0].map((_,i)=>a.map(r=>r[i]));
const mul=(a,b)=>a.map(r=>b[0].map((_,j)=>r.reduce((s,v,k)=>s+v*b[k][j],0)));
const mv=(a,v)=>a.map(r=>r.reduce((s,x,i)=>s+x*v[i],0));
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>Math.hypot(...a), unit=a=>a.map(x=>x/norm(a));
const labels=d=>d===2?['ux','uy','rz']:['ux','uy','uz','rx','ry','rz'];
function element(model,e){
 const d=model.dimension,p=d===2?3:6, a=model.nodes.find(n=>n.id===e.a),b=model.nodes.find(n=>n.id===e.b);
 if(!a||!b)throw Error(`Element ${e.id}: endpoint node does not exist.`);
 const delta=[b.x-a.x,b.y-a.y,d===3?b.z-a.z:0],L=norm(delta);
 if(!(L>0))throw Error(`Element ${e.id}: length must be positive.`);
 for(const key of d===2?['E','A','Iz']:['E','G','A','Iy','Iz','J'])if(!Number.isFinite(e[key])||e[key]<=0)throw Error(`Element ${e.id}: ${key} must be positive.`);
 const k=zeros(2*p,2*p),parts={axial:zeros(2*p,2*p),bendingZ:zeros(2*p,2*p)};
 function block(part,indices,values){indices.forEach((i,r)=>indices.forEach((j,c)=>{k[i][j]+=values[r][c];part[i][j]+=values[r][c];}));}
 const EA=e.E*e.A/L;
 block(parts.axial,[0,p],[[EA,-EA],[-EA,EA]]);
 const bending=(EI,sign)=>{const a=12*EI/L**3,b=sign*6*EI/L**2,c=4*EI/L,f=2*EI/L;return [[a,b,-a,b],[b,c,-b,f],[-a,-b,a,-b],[b,f,-b,c]];};
 block(parts.bendingZ,d===2?[1,2,4,5]:[1,5,7,11],bending(e.E*e.Iz,1));
 let R,T=zeros(2*p,2*p);
 if(d===2){const c=delta[0]/L,s=delta[1]/L;R=[[c,s,0],[-s,c,0],[0,0,1]];for(let h=0;h<2;h++)for(let i=0;i<3;i++)for(let j=0;j<3;j++)T[3*h+i][3*h+j]=R[i][j];}
 else {
  parts.bendingY=zeros(12,12);parts.torsion=zeros(12,12);
  block(parts.bendingY,[2,4,8,10],bending(e.E*e.Iy,-1));
  const GJ=e.G*e.J/L;block(parts.torsion,[3,9],[[GJ,-GJ],[-GJ,GJ]]);
  const x=unit(delta),ref=Math.abs(x[2])<0.9?[0,0,1]:[0,1,0],y0=unit(cross(ref,x)),z0=cross(x,y0),r=(e.roll||0)*Math.PI/180;
  const y=y0.map((v,i)=>v*Math.cos(r)+z0[i]*Math.sin(r)),z=z0.map((v,i)=>v*Math.cos(r)-y0[i]*Math.sin(r));R=[x,y,z];
  for(let h=0;h<4;h++)for(let i=0;i<3;i++)for(let j=0;j<3;j++)T[3*h+i][3*h+j]=R[i][j];
 }
 const ia=model.nodes.indexOf(a),ib=model.nodes.indexOf(b),map=[...Array(p).keys()].map(i=>ia*p+i).concat([...Array(p).keys()].map(i=>ib*p+i));
 return {e,L,k,parts,R,T,kg:mul(transpose(T),mul(k,T)),map,a,b};
}
function validate(m){
 if(![2,3].includes(m.dimension))throw Error('Dimension must be 2 or 3.');
 if(!Array.isArray(m.nodes)||!m.nodes.length||!Array.isArray(m.elements)||!m.elements.length)throw Error('Add at least two nodes and one element.');
 for(const key of ['nodes','elements']){const ids=m[key].map(x=>x.id);if(ids.some(x=>!Number.isInteger(x)||x<1)||new Set(ids).size!==ids.length)throw Error(`${key}: IDs must be unique positive integers.`);}
 const p=labels(m.dimension).length;
 for(const n of m.nodes){for(const key of m.dimension===2?['x','y']:['x','y','z'])if(!Number.isFinite(n[key]))throw Error(`Node ${n.id}: invalid ${key}.`);
 if(!Array.isArray(n.loads)||n.loads.length!==p||n.loads.some(x=>!Number.isFinite(x)))throw Error(`Node ${n.id}: provide ${p} finite loads.`);
 if(!Array.isArray(n.prescribed)||n.prescribed.length!==p||n.prescribed.some(x=>x!==null&&!Number.isFinite(x)))throw Error(`Node ${n.id}: restraints must be null (free) or a finite prescribed value.`);}
 for(const e of m.elements)if(e.roll!==undefined&&!Number.isFinite(e.roll))throw Error(`Element ${e.id}: roll must be finite.`);
}
function solve(m){
 validate(m);const p=labels(m.dimension).length,N=p*m.nodes.length,K=zeros(N,N),F=m.nodes.flatMap(n=>n.loads),u=Array(N).fill(0),free=[],fixed=[];
 const elements=m.elements.map(e=>element(m,e));
 for(const el of elements)el.map.forEach((i,r)=>el.map.forEach((j,c)=>K[i][j]+=el.kg[r][c]));
 m.nodes.forEach((n,i)=>n.prescribed.forEach((v,j)=>{if(v===null)free.push(i*p+j);else {fixed.push(i*p+j);u[i*p+j]=v;}}));
 const A=free.map(i=>free.map(j=>K[i][j])),rhs=free.map(i=>F[i]-fixed.reduce((s,j)=>s+K[i][j]*u[j],0));
 // Symmetric diagonal scaling controls translation/rotation unit disparity.
 const scales=A.map((r,i)=>{if(!(r[i]>0))throw Error(`Unstable structure: ${dofName(m,free[i])} has no stiffness.`);return Math.sqrt(r[i]);});
 const B=A.map((r,i)=>r.map((v,j)=>v/scales[i]/scales[j])),b=rhs.map((v,i)=>v/scales[i]),U=B.map(r=>r.slice()),g=b.slice(),steps=[];
 for(let k=0;k<free.length;k++){
  let pivot=k;for(let i=k+1;i<free.length;i++)if(Math.abs(U[i][k])>Math.abs(U[pivot][k]))pivot=i;
  if(Math.abs(U[pivot][k])<1e-12)throw Error('Unstable or numerically singular structure: add appropriate supports or check member connectivity and stiffness ratios.');
  if(pivot!==k){[U[k],U[pivot]]=[U[pivot],U[k]];[g[k],g[pivot]]=[g[pivot],g[k]];}
  const factors=[];for(let i=k+1;i<free.length;i++){const factor=U[i][k]/U[k][k];if(factor!==0)factors.push({row:i,factor});U[i][k]=0;for(let j=k+1;j<free.length;j++)U[i][j]-=factor*U[k][j];g[i]-=factor*g[k];}
  steps.push({column:k,pivotRow:pivot,pivot:U[k][k],factors});
 }
 const y=Array(free.length).fill(0),back=[];
 for(let i=free.length-1;i>=0;i--){const sum=U[i].reduce((s,v,j)=>j>i?s+v*y[j]:s,0);y[i]=(g[i]-sum)/U[i][i];back.push({row:i,rhs:g[i],sum,pivot:U[i][i],y:y[i]});}
 free.forEach((v,i)=>u[v]=y[i]/scales[i]);const Ku=mv(K,u),reactions=Ku.map((v,i)=>v-F[i]);
 if(u.some(v=>!Number.isFinite(v)))throw Error('Numerical overflow: check units and material properties.');
 elements.forEach(el=>{el.globalU=el.map.map(i=>u[i]);el.localU=mv(el.T,el.globalU);el.forces=mv(el.k,el.localU);});
 const residual=free.reduce((s,i)=>Math.max(s,Math.abs(reactions[i])),0),relativeResidual=residual/Math.max(1,...F.map(Math.abs),...Ku.map(Math.abs));
 return {K,F,u,free,fixed,A,rhs,scales,B,b,U,g,y,steps,back,reactions,elements,residual,relativeResidual};
}
function dofName(m,i){const p=labels(m.dimension).length;return `N${m.nodes[Math.floor(i/p)].id}.${labels(m.dimension)[i%p]}`;}
function curve(el,dimension,t){const v=el.localU,p=dimension===2?3:6,L=el.L,N1=1-3*t*t+2*t*t*t,N2=L*(t-2*t*t+t*t*t),N3=3*t*t-2*t*t*t,N4=L*(-t*t+t*t*t);let local;
 if(dimension===2)local=[v[0]*(1-t)+v[3]*t,N1*v[1]+N2*v[2]+N3*v[4]+N4*v[5],0];
 else local=[v[0]*(1-t)+v[6]*t,N1*v[1]+N2*v[5]+N3*v[7]+N4*v[11],N1*v[2]-N2*v[4]+N3*v[8]-N4*v[10]];
 return mv(transpose(el.R),local);
}
function node(id,x,y,z,p,loads,fix){return {id,x,y,z,loads:loads||Array(p).fill(0),prescribed:fix||Array(p).fill(null)};}
const examples={
 portal:()=>({name:'Textbook portal frame',dimension:2,lengthUnit:'in',forceUnit:'lb',nodes:[node(1,-60,0,0,3,null,[0,0,0]),node(2,-60,120,0,3,[10000,0,0]),node(3,60,120,0,3,[0,0,5000]),node(4,60,0,0,3,null,[0,0,0])],elements:[{id:1,a:1,b:2,E:30e6,A:10,Iz:200},{id:2,a:2,b:3,E:30e6,A:10,Iz:100},{id:3,a:3,b:4,E:30e6,A:10,Iz:200}]}),
 spatial:()=>({name:'3D cantilever • combined loading',dimension:3,lengthUnit:'m',forceUnit:'N',nodes:[node(1,0,0,0,6,null,[0,0,0,0,0,0]),node(2,3,0,0,6,[10000,-1000,-1500,600,0,0])],elements:[{id:1,a:1,b:2,E:210e9,G:80e9,A:0.01,Iy:8e-6,Iz:6e-6,J:1e-5,roll:0}]}),
 spaceFrame:()=>({name:'3D bent frame',dimension:3,lengthUnit:'m',forceUnit:'N',nodes:[node(1,0,0,0,6,null,[0,0,0,0,0,0]),node(2,0,0,3,6),node(3,3,0,3,6),node(4,3,2,3,6,[500,-800,-1200,200,0,0])],elements:[[1,2],[2,3],[3,4]].map(([a,b],i)=>({id:i+1,a,b,E:210e9,G:80e9,A:0.01,Iy:8e-6,Iz:6e-6,J:1e-5,roll:0}))})
};
const api={zeros,transpose,mul,mv,labels,element,solve,dofName,curve,examples};if(typeof module!=='undefined')module.exports=api;else root.FEM=api;
})(typeof window!=='undefined'?window:globalThis);
