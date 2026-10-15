/* repair.c - symmetry-free 2-opt repair of a near-diameter-4 4-regular graph.
 * usage: ./repair EDGEFILE SEED SECONDS BESTFILE
 * EDGEFILE: "a b" per line. BESTFILE: best edge list (resume input if exists).
 * defect = ordered pairs at distance > 4.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned long long u64;
typedef long long ll;
#define INF (1LL<<60)
#define MAXN 128
#define MAXE 256

static int n,M;
static int eu[MAXE],ev[MAXE];
static u64 abit[MAXN][2];
static u64 b2[MAXN][2];
static int defc[MAXN];
static u64 nmask0,nmask1;
static ll evals=0;
static u64 rs=1;
static u64 rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static int rndi(int k){return (int)(rnd()%(u64)k);}
static double now(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+1e-9*ts.tv_nsec;}

static inline int has(int a,int b){return (abit[a][b>>6]>>(b&63))&1ULL;}
static inline void setb(int a,int b){abit[a][b>>6]|=1ULL<<(b&63);}
static inline void clrb(int a,int b){abit[a][b>>6]&=~(1ULL<<(b&63));}

static ll evalg(void){
  evals++;
  for(int v=0;v<n;v++){
    u64 a0=abit[v][0],a1=abit[v][1],r0=a0,r1=a1,t;
    if(v<64)r0|=1ULL<<v;else r1|=1ULL<<(v-64);
    t=a0;while(t){int u=__builtin_ctzll(t);t&=t-1;r0|=abit[u][0];r1|=abit[u][1];}
    t=a1;while(t){int u=64+__builtin_ctzll(t);t&=t-1;r0|=abit[u][0];r1|=abit[u][1];}
    b2[v][0]=r0;b2[v][1]=r1;
  }
  ll def=0;
  for(int v=0;v<n;v++){
    u64 r0=0,r1=0,t;
    t=b2[v][0];while(t){int u=__builtin_ctzll(t);t&=t-1;r0|=b2[u][0];r1|=b2[u][1];}
    t=b2[v][1];while(t){int u=64+__builtin_ctzll(t);t&=t-1;r0|=b2[u][0];r1|=b2[u][1];}
    int cnt=n-__builtin_popcountll(r0&nmask0)-__builtin_popcountll(r1&nmask1);
    defc[v]=cnt;def+=cnt;
  }
  return def;
}

/* apply 2-opt on edges i,j with pattern p; returns 1 if applied (valid) */
static int apply2(int i,int j,int p){
  int a=eu[i],b=ev[i],c=eu[j],d=ev[j];
  if(p){int t=c;c=d;d=t;}
  /* new edges (a,c),(b,d) */
  if(a==c||b==d)return 0;
  if(a==b||c==d)return 0;
  if(a==d||b==c)return 0;
  if(has(a,c)||has(b,d))return 0;
  clrb(a,b);clrb(b,a);clrb(c,d);clrb(d,c);
  if(has(a,c)||has(b,d)){ /* (a,b) same as (c,d)? cannot: edges distinct */
    setb(a,b);setb(b,a);setb(c,d);setb(d,c);return 0;
  }
  setb(a,c);setb(c,a);setb(b,d);setb(d,b);
  eu[i]=a;ev[i]=c;eu[j]=b;ev[j]=d;
  return 1;
}
static void undo2(int i,int j,int a,int b,int c,int d){
  int na=eu[i],nb=ev[i],nc=eu[j],nd=ev[j];
  clrb(na,nb);clrb(nb,na);clrb(nc,nd);clrb(nd,nc);
  setb(a,b);setb(b,a);setb(c,d);setb(d,c);
  eu[i]=a;ev[i]=b;eu[j]=c;ev[j]=d;
}

/* all perfect matchings of 6 endpoints as index pairs into p[6] */
static const int m15[15][6]={
 {0,1,2,3,4,5},{0,1,2,4,3,5},{0,1,2,5,3,4},
 {0,2,1,3,4,5},{0,2,1,4,3,5},{0,2,1,5,3,4},
 {0,3,1,2,4,5},{0,3,1,4,2,5},{0,3,1,5,2,4},
 {0,4,1,2,3,5},{0,4,1,3,2,5},{0,4,1,5,2,3},
 {0,5,1,2,3,4},{0,5,1,3,2,4},{0,5,1,4,2,3}};
/* try all 3-opt rewirings of edges i,j,k; if best strictly < d, apply and return new d else d */
static ll try3(int i,int j,int k,ll d){
  int p[6];
  p[0]=eu[i];p[1]=ev[i];p[2]=eu[j];p[3]=ev[j];p[4]=eu[k];p[5]=ev[k];
  for(int x=0;x<6;x++)for(int y=x+1;y<6;y++)if(p[x]==p[y])return d;
  /* remove all three */
  clrb(p[0],p[1]);clrb(p[1],p[0]);clrb(p[2],p[3]);clrb(p[3],p[2]);clrb(p[4],p[5]);clrb(p[5],p[4]);
  int bm=-1;ll bd=d;
  for(int mi=1;mi<15;mi++){ /* mi=0 is identity */
    const int*q=m15[mi];
    int a1=p[q[0]],b1=p[q[1]],a2=p[q[2]],b2v=p[q[3]],a3=p[q[4]],b3=p[q[5]];
    if(has(a1,b1)||has(a2,b2v)||has(a3,b3))continue;
    setb(a1,b1);setb(b1,a1);
    if(has(a2,b2v)){clrb(a1,b1);clrb(b1,a1);continue;}
    setb(a2,b2v);setb(b2v,a2);
    if(has(a3,b3)){clrb(a1,b1);clrb(b1,a1);clrb(a2,b2v);clrb(b2v,a2);continue;}
    setb(a3,b3);setb(b3,a3);
    ll dd=evalg();
    if(dd<bd){bd=dd;bm=mi;}
    clrb(a1,b1);clrb(b1,a1);clrb(a2,b2v);clrb(b2v,a2);clrb(a3,b3);clrb(b3,a3);
  }
  if(bm>=0){
    const int*q=m15[bm];
    int a1=p[q[0]],b1=p[q[1]],a2=p[q[2]],b2v=p[q[3]],a3=p[q[4]],b3=p[q[5]];
    setb(a1,b1);setb(b1,a1);setb(a2,b2v);setb(b2v,a2);setb(a3,b3);setb(b3,a3);
    eu[i]=a1;ev[i]=b1;eu[j]=a2;ev[j]=b2v;eu[k]=a3;ev[k]=b3;
    return bd;
  }
  setb(p[0],p[1]);setb(p[1],p[0]);setb(p[2],p[3]);setb(p[3],p[2]);setb(p[4],p[5]);setb(p[5],p[4]);
  return d;
}
static int besteu[MAXE],bestev[MAXE];
static ll bstdef=INF;
static char bstf[300];
static void saveb(void){
  char tmp[320];snprintf(tmp,sizeof tmp,"%s.tmp",bstf);
  FILE*f=fopen(tmp,"w");if(!f)return;
  for(int i=0;i<M;i++)fprintf(f,"%d %d\n",besteu[i],bestev[i]);
  fclose(f);rename(tmp,bstf);
}

static int loadedges(const char*fn){
  FILE*f=fopen(fn,"r");if(!f)return 0;
  M=0;int a,b;n=0;
  while(fscanf(f,"%d %d",&a,&b)==2){eu[M]=a;ev[M]=b;M++;if(a>=n)n=a+1;if(b>=n)n=b+1;}
  fclose(f);
  memset(abit,0,sizeof abit);
  for(int i=0;i<M;i++){setb(eu[i],ev[i]);setb(ev[i],eu[i]);}
  return M>0;
}

int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: repair EDGEFILE SEED SECONDS BESTFILE\n");return 1;}
  rs=strtoull(argv[2],0,10)*2862933555777941757ULL+3037000493ULL;if(!rs)rs=1;
  double tlim=atof(argv[3]);
  snprintf(bstf,sizeof bstf,"%s",argv[4]);
  if(!loadedges(bstf)){
    if(!loadedges(argv[1])){fprintf(stderr,"no input\n");return 1;}
    printf("start from %s\n",argv[1]);
  }else printf("resume from %s\n",bstf);
  nmask0=n>=64?~0ULL:((1ULL<<n)-1);
  nmask1=n<=64?0ULL:(n>=128?~0ULL:((1ULL<<(n-64))-1));
  ll curdef=evalg();
  bstdef=curdef;memcpy(besteu,eu,sizeof(int)*M);memcpy(bestev,ev,sizeof(int)*M);
  printf("n=%d M=%d init def=%lld\n",n,M,curdef);fflush(stdout);
  if(argc>5&&argv[5][0]=='w'){ /* plateau walk on defect-0 graphs + insertion scans */
    double t0w=now(),tendw=t0w+tlim;
    ll wbest=INF;
    if(curdef!=0){printf("need defect-0 start\n");return 1;}
    while(now()<tendw){
      /* walk ~250 accepted 0-defect 2-opts */
      int acc=0;
      while(acc<250&&now()<tendw){
        int i=rndi(M),j=rndi(M);
        if(i==j)continue;
        int oa=eu[i],ob=ev[i],oc=eu[j],od=ev[j];
        if(oa==oc||oa==od||ob==oc||ob==od)continue;
        if(!apply2(i,j,rndi(2)))continue;
        if(evalg()==0)acc++;
        else undo2(i,j,oa,ob,oc,od);
      }
      /* insertion scan */
      int x=n;
      u64 sm0=nmask0,sm1=nmask1;
      n++;
      nmask0=n>=64?~0ULL:((1ULL<<n)-1);
      nmask1=n<=64?0ULL:(n>=128?~0ULL:((1ULL<<(n-64))-1));
      ll bd=INF;int bi=-1,bj=-1;
      for(int i=0;i<M;i++)for(int j=i+1;j<M;j++){
        int a=eu[i],b=ev[i],c=eu[j],d=ev[j];
        if(a==c||a==d||b==c||b==d)continue;
        clrb(a,b);clrb(b,a);clrb(c,d);clrb(d,c);
        setb(x,a);setb(a,x);setb(x,b);setb(b,x);setb(x,c);setb(c,x);setb(x,d);setb(d,x);
        ll dd=evalg();
        if(dd<bd){bd=dd;bi=i;bj=j;}
        clrb(x,a);clrb(a,x);clrb(x,b);clrb(b,x);clrb(x,c);clrb(c,x);clrb(x,d);clrb(d,x);
        setb(a,b);setb(b,a);setb(c,d);setb(d,c);
        if(bd==0)break;
        if(bd==0)break;
      }
      if(bd<wbest){
        wbest=bd;
        printf("t=%.1f walk insertion best=%lld\n",now()-t0w,wbest);fflush(stdout);
        FILE*f=fopen(bstf,"w");
        int a=eu[bi],b=ev[bi],c=eu[bj],d=ev[bj];
        for(int i=0;i<M;i++)if(i!=bi&&i!=bj)fprintf(f,"%d %d\n",eu[i],ev[i]);
        fprintf(f,"%d %d\n%d %d\n%d %d\n%d %d\n",x,a,x,b,x,c,x,d);
        fclose(f);
      }
      n--;nmask0=sm0;nmask1=sm1;
      if(wbest==0)break;
    }
    printf("RESULT walk n=%d best_insertion=%lld evals=%lld\n",n,wbest,evals);
    return 0;
  }
  if(argc>5&&argv[5][0]=='x'){ /* sampled general insertion: x->v1..v4, remove one incident edge each, rematch stubs */
    int x=n;n++;
    nmask0=n>=64?~0ULL:((1ULL<<n)-1);
    nmask1=n<=64?0ULL:(n>=128?~0ULL:((1ULL<<(n-64))-1));
    double t0x=now(),tendx=t0x+tlim;
    ll bd=INF;long long tried=0;
    while(now()<tendx){
      int v[4],ei[4],w[4];
      v[0]=rndi(x);v[1]=rndi(x);v[2]=rndi(x);v[3]=rndi(x);
      if(v[0]==v[1]||v[0]==v[2]||v[0]==v[3]||v[1]==v[2]||v[1]==v[3]||v[2]==v[3])continue;
      int ok=1;
      for(int a=0;a<4&&ok;a++){
        int cnt=0,pick=-1,seen=0;
        /* reservoir-pick a random incident edge */
        for(int e=0;e<M;e++)if(eu[e]==v[a]||ev[e]==v[a]){seen++;if(rndi(seen)==0)pick=e;}
        ei[a]=pick;
        w[a]=(eu[pick]==v[a])?ev[pick]:eu[pick];
        (void)cnt;
      }
      if(ei[0]==ei[1]||ei[0]==ei[2]||ei[0]==ei[3]||ei[1]==ei[2]||ei[1]==ei[3]||ei[2]==ei[3])continue;
      /* w's must not equal x-attachment conflicts: w can equal some v (fine unless creates dup) */
      tried++;
      /* remove 4 edges */
      for(int a=0;a<4;a++){clrb(eu[ei[a]],ev[ei[a]]);clrb(ev[ei[a]],eu[ei[a]]);}
      /* add x-v */
      int okx=1;
      for(int a=0;a<4;a++){if(has(x,v[a])){okx=0;break;}setb(x,v[a]);setb(v[a],x);}
      if(okx){
        static const int mt[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
        for(int s=0;s<3;s++){
          int a1=w[mt[s][0]],b1=w[mt[s][1]],a2=w[mt[s][2]],b2v=w[mt[s][3]];
          if(a1==b1||a2==b2v)continue;
          if(has(a1,b1))continue;
          setb(a1,b1);setb(b1,a1);
          if(has(a2,b2v)){clrb(a1,b1);clrb(b1,a1);continue;}
          setb(a2,b2v);setb(b2v,a2);
          ll dd=evalg();
          if(dd<bd){
            bd=dd;
            printf("t=%.1f xins def=%lld (tried %lld)\n",now()-t0x,bd,tried);fflush(stdout);
            FILE*f=fopen(bstf,"w");
            for(int e=0;e<M;e++)if(e!=ei[0]&&e!=ei[1]&&e!=ei[2]&&e!=ei[3])fprintf(f,"%d %d\n",eu[e],ev[e]);
            for(int a=0;a<4;a++)fprintf(f,"%d %d\n",x,v[a]);
            fprintf(f,"%d %d\n%d %d\n",a1,b1,a2,b2v);
            fclose(f);
          }
          clrb(a1,b1);clrb(b1,a1);clrb(a2,b2v);clrb(b2v,a2);
          if(bd==0)break;
        }
      }
      /* undo x-v and restore edges */
      for(int a=0;a<4;a++){clrb(x,v[a]);clrb(v[a],x);}
      for(int a=0;a<4;a++){setb(eu[ei[a]],ev[ei[a]]);setb(ev[ei[a]],eu[ei[a]]);}
      if(bd==0)break;
    }
    printf("RESULT xins n=%d best=%lld tried=%lld evals=%lld\n",n,bd,tried,evals);
    return 0;
  }
  if(argc>5&&argv[5][0]=='2'){ /* random-sampled adjacent-double-vertex insertion scan */
    int x=n,y=n+1;n+=2;
    nmask0=n>=64?~0ULL:((1ULL<<n)-1);
    nmask1=n<=64?0ULL:(n>=128?~0ULL:((1ULL<<(n-64))-1));
    double t02=now(),tend2=t02+tlim;
    ll bd=INF;long long tried=0;
    int bi=-1,bj=-1,bk=-1,bs=-1;
    static const int spl[10][6]={{0,1,2,3,4,5},{0,1,3,2,4,5},{0,1,4,2,3,5},{0,1,5,2,3,4},{0,2,3,1,4,5},{0,2,4,1,3,5},{0,2,5,1,3,4},{0,3,4,1,2,5},{0,3,5,1,2,4},{0,4,5,1,2,3}};
    while(now()<tend2){
      int i=rndi(M),j=rndi(M),k=rndi(M);
      if(i==j||j==k||i==k)continue;
      int p[6]={eu[i],ev[i],eu[j],ev[j],eu[k],ev[k]};
      int dup=0;
      for(int a2=0;a2<6;a2++)for(int b2=a2+1;b2<6;b2++)if(p[a2]==p[b2])dup=1;
      if(dup)continue;
      tried++;
      clrb(p[0],p[1]);clrb(p[1],p[0]);clrb(p[2],p[3]);clrb(p[3],p[2]);clrb(p[4],p[5]);clrb(p[5],p[4]);
      setb(x,y);setb(y,x);
      for(int s=0;s<10;s++){
        const int*q=spl[s];
        int xa=p[q[0]],xb=p[q[1]],xc=p[q[2]],ya=p[q[3]],yb=p[q[4]],yc=p[q[5]];
        setb(x,xa);setb(xa,x);setb(x,xb);setb(xb,x);setb(x,xc);setb(xc,x);
        setb(y,ya);setb(ya,y);setb(y,yb);setb(yb,y);setb(y,yc);setb(yc,y);
        ll dd=evalg();
        if(dd<bd){
          bd=dd;bi=i;bj=j;bk=k;bs=s;
          printf("t=%.1f dbl-ins def=%lld (tried %lld)\n",now()-t02,bd,tried);fflush(stdout);
          FILE*f=fopen(bstf,"w");
          for(int e=0;e<M;e++)if(e!=i&&e!=j&&e!=k)fprintf(f,"%d %d\n",eu[e],ev[e]);
          fprintf(f,"%d %d\n",x,y);
          fprintf(f,"%d %d\n%d %d\n%d %d\n",x,xa,x,xb,x,xc);
          fprintf(f,"%d %d\n%d %d\n%d %d\n",y,ya,y,yb,y,yc);
          fclose(f);
        }
        clrb(x,xa);clrb(xa,x);clrb(x,xb);clrb(xb,x);clrb(x,xc);clrb(xc,x);
        clrb(y,ya);clrb(ya,y);clrb(y,yb);clrb(yb,y);clrb(y,yc);clrb(yc,y);
        if(bd==0)break;
      }
      clrb(x,y);clrb(y,x);
      setb(p[0],p[1]);setb(p[1],p[0]);setb(p[2],p[3]);setb(p[3],p[2]);setb(p[4],p[5]);setb(p[5],p[4]);
      if(bd==0)break;
    }
    printf("RESULT dbl n=%d best=%lld tried=%lld evals=%lld\n",n,bd,tried,evals);
    return 0;
  }
  if(argc>5&&argv[5][0]=='g'){ /* exhaustive single-vertex insertion scan */
    int x=n;n++;
    nmask0=n>=64?~0ULL:((1ULL<<n)-1);
    nmask1=n<=64?0ULL:(n>=128?~0ULL:((1ULL<<(n-64))-1));
    ll bd=INF;int bi=-1,bj=-1;
    for(int i=0;i<M;i++)for(int j=i+1;j<M;j++){
      int a=eu[i],b=ev[i],c=eu[j],d=ev[j];
      if(a==c||a==d||b==c||b==d)continue;
      clrb(a,b);clrb(b,a);clrb(c,d);clrb(d,c);
      setb(x,a);setb(a,x);setb(x,b);setb(b,x);setb(x,c);setb(c,x);setb(x,d);setb(d,x);
      ll dd=evalg();
      if(dd<bd){bd=dd;bi=i;bj=j;printf("ins (%d,%d)(%d,%d) def=%lld\n",a,b,c,d,dd);fflush(stdout);}
      clrb(x,a);clrb(a,x);clrb(x,b);clrb(b,x);clrb(x,c);clrb(c,x);clrb(x,d);clrb(d,x);
      setb(a,b);setb(b,a);setb(c,d);setb(d,c);
    }
    printf("BEST insertion def=%lld edges %d,%d\n",bd,bi,bj);
    /* write best insertion */
    {
      int a=eu[bi],b=ev[bi],c=eu[bj],d=ev[bj];
      FILE*f=fopen(bstf,"w");
      for(int i=0;i<M;i++)if(i!=bi&&i!=bj)fprintf(f,"%d %d\n",eu[i],ev[i]);
      fprintf(f,"%d %d\n%d %d\n%d %d\n%d %d\n",x,a,x,b,x,c,x,d);
      fclose(f);
    }
    return 0;
  }
  if(argc>5&&argv[5][0]=='a'){ /* simulated annealing */
    double T0=argc>6?atof(argv[6]):3.0;
    double T=T0;
    double t0a=now(),tenda=t0a+tlim;
    ll steps=0;
    while(bstdef>0&&now()<tenda){
      steps++;
      T*=0.9999985;
      if(T<0.15)T=T0*(0.3+0.07*rndi(10));
      int i=rndi(M),j=rndi(M);
      if(i==j)continue;
      int oa=eu[i],ob=ev[i],oc=eu[j],od=ev[j];
      if(oa==oc||oa==od||ob==oc||ob==od)continue;
      int p=rndi(2);
      if(!apply2(i,j,p))continue;
      ll d=evalg();
      ll dl=d-curdef;
      if(dl<=0||rnd()%1000000<(u64)(1000000.0*exp(-(double)dl/T))){
        curdef=d;
        if(curdef<bstdef){
          bstdef=curdef;memcpy(besteu,eu,sizeof(int)*M);memcpy(bestev,ev,sizeof(int)*M);
          printf("t=%.1f best=%lld T=%.3f steps=%lld\n",now()-t0a,bstdef,T,steps);fflush(stdout);
          saveb();
        }
      }else undo2(i,j,oa,ob,oc,od);
    }
    saveb();
    printf("RESULT n=%d best=%lld evals=%lld (anneal)\n",n,bstdef,evals);
    return 0;
  }
  if(tlim<0){ /* exhaustive targeted 3-opt sweep */
    double tend2=now()-tlim;
    int changed=1;
    while(changed&&now()<tend2&&curdef>0){
      changed=0;
      evalg();
      int tset[40];int nt=0;
      for(int i=0;i<M&&nt<40;i++)if(defc[eu[i]]>0||defc[ev[i]]>0)tset[nt++]=i;
      printf("sweep: def=%lld targeted_edges=%d\n",curdef,nt);fflush(stdout);
      for(int a=0;a<nt&&!changed;a++){
        for(int j=0;j<M&&!changed;j++){
          if(j==tset[a])continue;
          for(int k=j+1;k<M;k++){
            if(k==tset[a])continue;
            ll nd3=try3(tset[a],j,k,curdef);
            if(nd3<curdef){
              curdef=nd3;changed=1;
              printf("3opt -> %lld evals=%lld\n",curdef,evals);fflush(stdout);
              if(curdef<bstdef){bstdef=curdef;memcpy(besteu,eu,sizeof(int)*M);memcpy(bestev,ev,sizeof(int)*M);saveb();}
              break;
            }
            if(((k)&511)==0&&now()>tend2){changed=0;a=nt;j=M;break;}
          }
        }
      }
    }
    printf("RESULT n=%d best=%lld evals=%lld (3opt-sweep)\n",n,bstdef,evals);
    return 0;
  }
  double t0=now(),tend=t0+tlim;
  int stall=0;
  while(bstdef>0&&now()<tend){
    /* targeted step: pick deficient vertex, incident edge, scan all partners */
    ll tot=0;for(int v=0;v<n;v++)tot+=defc[v];
    int u;
    if(tot<=0)break;
    if(rndi(100)<15)u=rndi(n);
    else{ll p=(ll)(rnd()%(u64)tot);int v=0;while(p>=defc[v]){p-=defc[v];v++;}u=v;}
    int inc[8],ni=0;
    for(int i=0;i<M;i++)if(eu[i]==u||ev[i]==u){inc[ni++]=i;if(ni>=8)break;}
    if(!ni){stall++;continue;}
    int e1=inc[rndi(ni)];
    int ba=eu[e1],bb=ev[e1];
    ll bd=curdef;int bj=-1,bp=0;
    int eqj[MAXE],eqp[MAXE],neq=0;
    for(int j=0;j<M;j++){
      if(j==e1)continue;
      int ca=eu[j],cb=ev[j];
      if(ca==ba||ca==bb||cb==ba||cb==bb)continue;
      for(int p=0;p<2;p++){
        int oa=eu[e1],ob=ev[e1],oc=eu[j],od=ev[j];
        if(!apply2(e1,j,p))continue;
        ll d=evalg();
        if(d<bd){bd=d;bj=j;bp=p;}
        else if(d==curdef&&neq<MAXE){eqj[neq]=j;eqp[neq]=p;neq++;}
        undo2(e1,j,oa,ob,oc,od);
      }
    }
    if(bj>=0){
      apply2(e1,bj,bp);
      curdef=evalg();
      if(curdef<bstdef){
        bstdef=curdef;memcpy(besteu,eu,sizeof(int)*M);memcpy(bestev,ev,sizeof(int)*M);
        printf("t=%.1f best=%lld evals=%lld\n",now()-t0,bstdef,evals);fflush(stdout);
        saveb();
      }
      stall=0;
    }else{
      if(neq>0&&rndi(100)<40){
        int k=rndi(neq);
        apply2(e1,eqj[k],eqp[k]);
        curdef=evalg();
      }
      stall++;
    }
    if(stall==61&&curdef<=24&&curdef<=bstdef+4){
      /* exhaustive targeted 3-opt: edges incident to deficient vertices vs all pairs */
      evalg();
      int tset[16];int nt=0;
      for(int i=0;i<M&&nt<16;i++)if(defc[eu[i]]>0||defc[ev[i]]>0)tset[nt++]=i;
      ll d3=curdef;int done=0;
      for(int a=0;a<nt&&!done;a++){
        for(int j=0;j<M&&!done;j++){
          if(j==tset[a])continue;
          for(int k=j+1;k<M;k++){
            if(k==tset[a])continue;
            ll nd3=try3(tset[a],j,k,d3);
            if(nd3<d3){
              d3=nd3;
              if(now()>tend){done=1;break;}
            }
            if(((j*M+k)&1023)==0&&now()>tend){done=1;break;}
          }
        }
      }
      curdef=evalg();
      if(curdef<bstdef){
        bstdef=curdef;memcpy(besteu,eu,sizeof(int)*M);memcpy(bestev,ev,sizeof(int)*M);
        printf("t=%.1f best=%lld evals=%lld (3opt)\n",now()-t0,bstdef,evals);fflush(stdout);
        saveb();stall=0;continue;
      }
    }
    if(stall>60){
      /* kick: restart from best (75%) + random 2-opts */
      if(rndi(100)<75){
        memcpy(eu,besteu,sizeof(int)*M);memcpy(ev,bestev,sizeof(int)*M);
        memset(abit,0,sizeof abit);
        for(int i=0;i<M;i++){setb(eu[i],ev[i]);setb(ev[i],eu[i]);}
      }
      int k=1+rndi(3);
      for(int s=0;s<k;s++){
        for(int att=0;att<50;att++){
          int i=rndi(M),j=rndi(M);
          if(i==j)continue;
          int ca=eu[j],cb=ev[j],oa=eu[i],ob=ev[i];
          if(ca==oa||ca==ob||cb==oa||cb==ob)continue;
          if(apply2(i,j,rndi(2)))break;
        }
      }
      curdef=evalg();
      stall=0;
    }
  }
  saveb();
  printf("RESULT n=%d best=%lld evals=%lld\n",n,bstdef,evals);
  return 0;
}
