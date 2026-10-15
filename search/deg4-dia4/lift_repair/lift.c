/* lift.c - voltage-graph lift search for 4-regular diameter<=4 graphs.
 * Supports arbitrary finite groups via multiplication tables (cyclic,
 * dihedral, dicyclic, Frobenius p:q:r) plus fixed vertices (FIXO classes,
 * requires |G|=4) for mixed cycle-type constructions.
 *
 * usage: ./lift GROUP B F SEED SECONDS CKPTFILE
 *   GROUP: cyc:K | dih:K (order 2K) | dic:K (order 4K) | frob:p:q:r
 *   SECONDS=0 -> load ckpt, dump edge list to CKPTFILE.edges, exit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;
typedef long long ll;
#define INF (1LL<<60)
#define MAXM 60
#define MAXN 128
#define MAXB 34
#define MAXC 64

enum{LINK,LOOP,SEMI,FIXO};
typedef struct{int t,u,v,g;}Cls;

static int m,B,F,N;
static int mul[MAXM][MAXM],ginv[MAXM];
static int invl[MAXM],ninv,loopl[MAXM],nlp;
static u64 nmask0,nmask1;
static Cls cur[MAXC],bst[MAXC];
static int ncur,nbst;
static ll bstdef=INF,evals=0;
static int repdef[MAXB+8];
static char gname[64],ckptf[300];
static u64 rs=1;

static u64 rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static int rndi(int k){return (int)(rnd()%(u64)k);}
static double now(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+1e-9*ts.tv_nsec;}

static void build_group(const char*s){
  int i,j;
  if(!strncmp(s,"cyc:",4)){int k=atoi(s+4);m=k;
    for(i=0;i<k;i++)for(j=0;j<k;j++)mul[i][j]=(i+j)%k;
  }else if(!strncmp(s,"dih:",4)){int k=atoi(s+4);m=2*k;
    for(int t1=0;t1<2;t1++)for(int a1=0;a1<k;a1++)for(int t2=0;t2<2;t2++)for(int a2=0;a2<k;a2++){
      int t=t1^t2;int a=(t2?(k-a1)%k:a1);a=(a+a2)%k;
      mul[t1*k+a1][t2*k+a2]=t*k+a;
    }
  }else if(!strncmp(s,"dic:",4)){int k=atoi(s+4);m=4*k;int K=2*k;
    for(int t1=0;t1<2;t1++)for(int a1=0;a1<K;a1++)for(int t2=0;t2<2;t2++)for(int a2=0;a2<K;a2++){
      int t,a;
      if(!t1){t=t2;a=(a1+a2)%K;}
      else{int d=(a1-a2+K)%K;if(!t2){t=1;a=d;}else{t=0;a=(d+k)%K;}}
      mul[t1*K+a1][t2*K+a2]=t*K+a;
    }
  }else if(!strncmp(s,"frob:",5)){int p,q,r;
    if(sscanf(s+5,"%d:%d:%d",&p,&q,&r)!=3){fprintf(stderr,"bad frob\n");exit(1);}
    m=p*q;
    int pw[16];pw[0]=1;for(i=1;i<q;i++)pw[i]=pw[i-1]*r%p;
    if(pw[q-1]*r%p!=1){fprintf(stderr,"frob: r order wrong\n");exit(1);}
    for(int y1=0;y1<q;y1++)for(int x1=0;x1<p;x1++)for(int y2=0;y2<q;y2++)for(int x2=0;x2<p;x2++)
      mul[y1*p+x1][y2*p+x2]=((y1+y2)%q)*p+(x1+pw[y1]*x2)%p;
  }else if(!strncmp(s,"prod:",5)){int a,b;
    if(sscanf(s+5,"%d:%d",&a,&b)!=2){fprintf(stderr,"bad prod\n");exit(1);}
    m=a*b;
    for(int x1=0;x1<a;x1++)for(int y1=0;y1<b;y1++)for(int x2=0;x2<a;x2++)for(int y2=0;y2<b;y2++)
      mul[y1*a+x1][y2*a+x2]=((y1+y2)%b)*a+(x1+x2)%a;
  }else if(!strcmp(s,"e8")){m=8;
    for(i=0;i<8;i++)for(j=0;j<8;j++)mul[i][j]=i^j;
  }else if(!strcmp(s,"v4")){m=4;
    for(i=0;i<4;i++)for(j=0;j<4;j++)mul[i][j]=i^j;
  }else if(!strcmp(s,"a4")){
    int perms[12][4],np=0;
    for(int a=0;a<4;a++)for(int b=0;b<4;b++){if(b==a)continue;
      for(int c=0;c<4;c++){if(c==a||c==b)continue;int d=6-a-b-c;
        int q[4];q[0]=a;q[1]=b;q[2]=c;q[3]=d;
        int iv=0;for(int x=0;x<4;x++)for(int y=x+1;y<4;y++)if(q[x]>q[y])iv++;
        if(iv%2==0){memcpy(perms[np],q,sizeof q);np++;}
      }}
    m=12;
    for(i=0;i<12;i++)for(j=0;j<12;j++){
      int comp[4];for(int x=0;x<4;x++)comp[x]=perms[i][perms[j][x]];
      for(int k2=0;k2<12;k2++)if(!memcmp(comp,perms[k2],sizeof comp)){mul[i][j]=k2;break;}
    }
  }else{fprintf(stderr,"bad group %s\n",s);exit(1);}
  if(m>MAXM){fprintf(stderr,"m too big\n");exit(1);}
  for(i=0;i<m;i++)if(mul[0][i]!=i||mul[i][0]!=i){fprintf(stderr,"identity broken\n");exit(1);}
  /* associativity spot check */
  for(int c=0;c<200;c++){int a=rndi(m),b=rndi(m),d=rndi(m);
    if(mul[mul[a][b]][d]!=mul[a][mul[b][d]]){fprintf(stderr,"assoc broken\n");exit(1);}}
  ninv=nlp=0;
  for(i=0;i<m;i++){ginv[i]=-1;for(j=0;j<m;j++)if(mul[i][j]==0){ginv[i]=j;break;}
    if(ginv[i]<0){fprintf(stderr,"no inverse\n");exit(1);}}
  for(i=1;i<m;i++){if(mul[i][i]==0)invl[ninv++]=i;else if(i<=ginv[i])loopl[nlp++]=i;}
}

static u64 adjb[MAXN][2],b2[MAXN][2];
static inline int ae(int a,int b){
  if(a==b)return 1;
  if((adjb[a][b>>6]>>(b&63))&1ULL)return 1;
  adjb[a][b>>6]|=1ULL<<(b&63);adjb[b][a>>6]|=1ULL<<(a&63);return 0;
}
static ll evalcfg(Cls*C,int nc){
  evals++;
  memset(adjb,0,sizeof(u64)*2*(size_t)N);
  int bad=0;
  for(int i=0;i<nc;i++){
    Cls*c=&C[i];
    if(c->t==LINK){
      if(c->u==c->v)return INF;
      int base=c->v*m,ub=c->u*m;
      for(int x=0;x<m;x++)bad|=ae(ub+x,base+mul[x][c->g]);
    }else if(c->t==LOOP){
      if(c->g==0||mul[c->g][c->g]==0)return INF;
      int base=c->u*m;
      for(int x=0;x<m;x++)bad|=ae(base+x,base+mul[x][c->g]);
    }else if(c->t==SEMI){
      if(c->g==0||mul[c->g][c->g]!=0)return INF;
      int base=c->u*m;
      for(int x=0;x<m;x++){int y=mul[x][c->g];if(x<y)bad|=ae(base+x,base+y);}
    }else{
      int fv=B*m+c->u,base=c->v*m;
      for(int x=0;x<m;x++)bad|=ae(fv,base+x);
    }
    if(bad)return INF;
  }
  for(int v=0;v<N;v++){
    int d=__builtin_popcountll(adjb[v][0])+__builtin_popcountll(adjb[v][1]);
    if(d!=4)return INF;
  }
  for(int v=0;v<N;v++){
    u64 a0=adjb[v][0],a1=adjb[v][1],r0=a0,r1=a1,t;
    if(v<64)r0|=1ULL<<v;else r1|=1ULL<<(v-64);
    t=a0;while(t){int u=__builtin_ctzll(t);t&=t-1;r0|=adjb[u][0];r1|=adjb[u][1];}
    t=a1;while(t){int u=64+__builtin_ctzll(t);t&=t-1;r0|=adjb[u][0];r1|=adjb[u][1];}
    b2[v][0]=r0;b2[v][1]=r1;
  }
  ll def=0;
  for(int ri=0;ri<B+F;ri++){
    int rep=ri<B?ri*m:B*m+(ri-B);
    int mult=ri<B?m:1;
    u64 r0=0,r1=0,t;
    t=b2[rep][0];while(t){int u=__builtin_ctzll(t);t&=t-1;r0|=b2[u][0];r1|=b2[u][1];}
    t=b2[rep][1];while(t){int u=64+__builtin_ctzll(t);t&=t-1;r0|=b2[u][0];r1|=b2[u][1];}
    int cnt=N-__builtin_popcountll(r0&nmask0)-__builtin_popcountll(r1&nmask1);
    repdef[ri]=cnt;def+=(ll)mult*cnt;
  }
  return def;
}

static int allow_side=0;
static ll scan_class(Cls*C,int nc,int ci,ll d0){
  Cls*c=&C[ci];
  if(c->t==FIXO)return d0;
  int bg=c->g;ll bd=d0;
  const int*lst;int nl;
  static int all[MAXM];
  int eq[MAXM],neq=0;
  if(c->t==LINK){for(int g=0;g<m;g++)all[g]=g;lst=all;nl=m;}
  else if(c->t==LOOP){lst=loopl;nl=nlp;}
  else{lst=invl;nl=ninv;}
  int og=c->g;
  for(int k=0;k<nl;k++){
    int g=lst[k];if(g==og)continue;
    c->g=g;ll d=evalcfg(C,nc);
    if(d<bd){bd=d;bg=g;}
    else if(d==d0&&allow_side)eq[neq++]=g;
  }
  if(bd>=d0&&neq>0&&allow_side&&rndi(100)<35)bg=eq[rndi(neq)];
  c->g=bg;
  return bd;
}

static ll conv_polish(ll d,int*imp){
  Cls bak[MAXC];int nbak;
  if(ninv>0)for(int i=0;i<ncur&&ncur<MAXC;i++)if(cur[i].t==LINK){
    memcpy(bak,cur,sizeof(Cls)*ncur);nbak=ncur;
    int u=cur[i].u,v=cur[i].v;
    cur[i].t=SEMI;cur[i].u=u;cur[i].v=0;cur[i].g=invl[0];
    cur[ncur].t=SEMI;cur[ncur].u=v;cur[ncur].v=0;cur[ncur].g=invl[0];ncur++;
    ll d1=evalcfg(cur,ncur);
    d1=scan_class(cur,ncur,i,d1);
    d1=scan_class(cur,ncur,ncur-1,d1);
    d1=scan_class(cur,ncur,i,d1);
    if(d1<d){d=d1;*imp=1;}
    else{memcpy(cur,bak,sizeof(Cls)*nbak);ncur=nbak;}
  }
  for(int i=0;i<ncur;i++)if(cur[i].t==SEMI)for(int j=i+1;j<ncur;j++){
    if(cur[i].t!=SEMI)break;
    if(cur[j].t!=SEMI)continue;
    int u=cur[i].u,v=cur[j].u;
    if(u==v&&nlp==0)continue;
    memcpy(bak,cur,sizeof(Cls)*ncur);nbak=ncur;
    if(u!=v){cur[i].t=LINK;cur[i].u=u;cur[i].v=v;cur[i].g=0;}
    else{cur[i].t=LOOP;cur[i].g=loopl[0];}
    cur[j]=cur[ncur-1];ncur--;
    ll d1=evalcfg(cur,ncur);
    d1=scan_class(cur,ncur,i,d1);
    if(d1<d){d=d1;*imp=1;break;}
    else{memcpy(cur,bak,sizeof(Cls)*nbak);ncur=nbak;}
  }
  if(ninv>=2)for(int i=0;i<ncur&&ncur<MAXC;i++)if(cur[i].t==LOOP){
    memcpy(bak,cur,sizeof(Cls)*ncur);nbak=ncur;
    int u=cur[i].u;
    cur[i].t=SEMI;cur[i].g=invl[0];
    cur[ncur].t=SEMI;cur[ncur].u=u;cur[ncur].v=0;cur[ncur].g=invl[1];ncur++;
    ll d1=evalcfg(cur,ncur);
    d1=scan_class(cur,ncur,i,d1);
    d1=scan_class(cur,ncur,ncur-1,d1);
    d1=scan_class(cur,ncur,i,d1);
    if(d1<d){d=d1;*imp=1;}
    else{memcpy(cur,bak,sizeof(Cls)*nbak);ncur=nbak;}
  }
  return d;
}

static ll polish(ll d,int deep){
  int imp=1;
  while(imp){
    imp=0;
    for(int ci=0;ci<ncur;ci++){
      ll nd=scan_class(cur,ncur,ci,d);
      if(nd<d){d=nd;imp=1;}
    }
    if(deep){
      for(int i=0;i<ncur;i++)if(cur[i].t==LINK)for(int j=i+1;j<ncur;j++)if(cur[j].t==LINK){
        for(int pat=0;pat<2;pat++){
          Cls s1=cur[i],s2=cur[j];
          if(pat==0){int t=cur[i].v;cur[i].v=cur[j].u;cur[j].u=t;}
          else{int t=cur[i].v;cur[i].v=cur[j].v;cur[j].v=t;}
          if(cur[i].u==cur[i].v||cur[j].u==cur[j].v){cur[i]=s1;cur[j]=s2;continue;}
          ll d1=evalcfg(cur,ncur);
          d1=scan_class(cur,ncur,i,d1);
          d1=scan_class(cur,ncur,j,d1);
          if(d1<d){d=d1;imp=1;}
          else{cur[i]=s1;cur[j]=s2;}
        }
      }
      d=conv_polish(d,&imp);
      if(!imp){
        /* coupled two-class voltage scan */
        for(int i=0;i<ncur;i++){
          if(cur[i].t==FIXO)continue;
          const int*li;int nli;
          static int alli[MAXM];
          if(cur[i].t==LINK){for(int g=0;g<m;g++)alli[g]=g;li=alli;nli=m;}
          else if(cur[i].t==LOOP){li=loopl;nli=nlp;}
          else{li=invl;nli=ninv;}
          for(int j=i+1;j<ncur;j++){
            if(cur[j].t==FIXO)continue;
            const int*lj;int nlj;
            static int allj[MAXM];
            if(cur[j].t==LINK){for(int g=0;g<m;g++)allj[g]=g;lj=allj;nlj=m;}
            else if(cur[j].t==LOOP){lj=loopl;nlj=nlp;}
            else{lj=invl;nlj=ninv;}
            if((ll)nli*nlj>2500)continue;
            int ogi=cur[i].g,ogj=cur[j].g,bgi=ogi,bgj=ogj;ll bd=d;
            for(int a=0;a<nli;a++){
              cur[i].g=li[a];
              for(int b=0;b<nlj;b++){
                if(li[a]==ogi&&lj[b]==ogj)continue;
                cur[j].g=lj[b];
                ll dd=evalcfg(cur,ncur);
                if(dd<bd){bd=dd;bgi=li[a];bgj=lj[b];}
              }
            }
            cur[i].g=bgi;cur[j].g=bgj;
            if(bd<d){d=bd;imp=1;}
          }
        }
      }
    }
  }
  return d;
}

static void rand_struct(void){
  for(int tries=0;tries<20;tries++){
    int op=rndi(F>0?5:4);
    if(op==0){
      int li[MAXC],nl=0;
      for(int i=0;i<ncur;i++)if(cur[i].t==LINK)li[nl++]=i;
      if(nl<2)continue;
      int i=li[rndi(nl)],j=li[rndi(nl)];
      if(i==j)continue;
      Cls s1=cur[i],s2=cur[j];
      if(rndi(2)){int t=cur[i].v;cur[i].v=cur[j].u;cur[j].u=t;}
      else{int t=cur[i].v;cur[i].v=cur[j].v;cur[j].v=t;}
      cur[i].g=rndi(m);cur[j].g=rndi(m);
      if(cur[i].u==cur[i].v||cur[j].u==cur[j].v){cur[i]=s1;cur[j]=s2;continue;}
      return;
    }else if(op==1&&ninv>0&&ncur<MAXC){
      int li[MAXC],nl=0;
      for(int i=0;i<ncur;i++)if(cur[i].t==LINK)li[nl++]=i;
      if(!nl)continue;
      int i=li[rndi(nl)];
      int u=cur[i].u,v=cur[i].v;
      cur[i].t=SEMI;cur[i].u=u;cur[i].v=0;cur[i].g=invl[rndi(ninv)];
      cur[ncur].t=SEMI;cur[ncur].u=v;cur[ncur].v=0;cur[ncur].g=invl[rndi(ninv)];ncur++;
      return;
    }else if(op==2){
      int si[MAXC],ns=0;
      for(int i=0;i<ncur;i++)if(cur[i].t==SEMI)si[ns++]=i;
      if(ns<2)continue;
      int i=si[rndi(ns)],j=si[rndi(ns)];
      if(i==j||cur[i].u==cur[j].u)continue;
      cur[i].t=LINK;cur[i].v=cur[j].u;cur[i].g=rndi(m);
      cur[j]=cur[ncur-1];ncur--;
      return;
    }else if(op==3){
      if(rndi(2)&&ninv>=2&&ncur<MAXC){
        int li[MAXC],nl=0;
        for(int i=0;i<ncur;i++)if(cur[i].t==LOOP)li[nl++]=i;
        if(!nl)continue;
        int i=li[rndi(nl)];int u=cur[i].u;
        cur[i].t=SEMI;cur[i].g=invl[rndi(ninv)];
        cur[ncur].t=SEMI;cur[ncur].u=u;cur[ncur].v=0;cur[ncur].g=invl[rndi(ninv)];ncur++;
        return;
      }else if(nlp>0){
        int si[MAXC],ns=0;
        for(int i=0;i<ncur;i++)if(cur[i].t==SEMI)si[ns++]=i;
        int fi=-1,fj=-1;
        for(int a=0;a<ns&&fi<0;a++)for(int b=a+1;b<ns;b++)if(cur[si[a]].u==cur[si[b]].u){fi=si[a];fj=si[b];break;}
        if(fi<0)continue;
        int u=cur[fi].u;
        cur[fi].t=LOOP;cur[fi].u=u;cur[fi].v=0;cur[fi].g=loopl[rndi(nlp)];
        cur[fj]=cur[ncur-1];ncur--;
        return;
      }else continue;
    }else if(op==4&&F>0){
      int fi2[16],nf=0;
      for(int i=0;i<ncur;i++)if(cur[i].t==FIXO)fi2[nf++]=i;
      if(!nf)continue;
      if(nf>=2&&rndi(2)){
        int i=fi2[rndi(nf)],j=fi2[rndi(nf)];
        if(i==j||cur[i].v==cur[j].v)continue;
        int t=cur[i].v;cur[i].v=cur[j].v;cur[j].v=t;return;
      }
      int i=fi2[rndi(nf)];int uo=cur[i].v;
      int j=rndi(ncur);Cls*c=&cur[j];
      if(c->t==LINK){
        int side=rndi(2);int w=side?c->v:c->u;int o=side?c->u:c->v;
        if(w==uo||o==uo)continue;
        if(side)c->v=uo;else c->u=uo;
        cur[i].v=w;return;
      }else if(c->t==SEMI){
        int w=c->u;if(w==uo)continue;
        c->u=uo;cur[i].v=w;return;
      }else continue;
    }
  }
}

static int init_struct(void){
  ncur=0;
  int defc[MAXB];
  for(int u=0;u<B;u++)defc[u]=4;
  for(int j=0;j<F;j++){
    if(ncur>=MAXC)return 0;
    int tgt=j%B;
    Cls c={FIXO,j,tgt,0};cur[ncur++]=c;
    defc[tgt]--;if(defc[tgt]<0)return 0;
  }
  for(int u=0;u+1<B;u++){
    if(defc[u]<1||defc[u+1]<1)return 0;
    if(ncur>=MAXC)return 0;
    Cls c={LINK,u,u+1,rndi(m)};cur[ncur++]=c;
    defc[u]--;defc[u+1]--;
  }
  while(1){
    int u=-1,mx=0;
    for(int v2=0;v2<B;v2++)if(defc[v2]>mx){mx=defc[v2];u=v2;}
    if(u<0)break;
    int cand[MAXB],ncand=0;
    for(int v2=0;v2<B;v2++)if(v2!=u&&defc[v2]>0)cand[ncand++]=v2;
    int choice=-1;
    if(ncand>0)choice=0;
    if(defc[u]>=2&&nlp>0&&rndi(100)<25)choice=1;
    if(ninv>0&&rndi(100)<12)choice=2;
    if(choice==0&&ncand==0)choice=-1;
    if(choice==1&&!(defc[u]>=2&&nlp>0))choice=-1;
    if(choice==2&&ninv==0)choice=-1;
    if(choice<0){
      if(ncand>0)choice=0;
      else if(defc[u]>=2&&nlp>0)choice=1;
      else if(ninv>0)choice=2;
      else return 0;
    }
    if(ncur>=MAXC)return 0;
    if(choice==0){
      int v2=cand[rndi(ncand)];
      Cls c={LINK,u,v2,rndi(m)};cur[ncur++]=c;
      defc[u]--;defc[v2]--;
    }else if(choice==1){
      Cls c={LOOP,u,0,loopl[rndi(nlp)]};cur[ncur++]=c;
      defc[u]-=2;
    }else{
      Cls c={SEMI,u,0,invl[rndi(ninv)]};cur[ncur++]=c;
      defc[u]--;
    }
  }
  return 1;
}

static void save_ckpt(void){
  char tmp[320];snprintf(tmp,sizeof tmp,"%s.tmp",ckptf);
  FILE*f=fopen(tmp,"w");if(!f)return;
  fprintf(f,"%lld %lld %d %d %d %d\n",bstdef,evals,m,B,F,nbst);
  for(int i=0;i<nbst;i++)fprintf(f,"%d %d %d %d\n",bst[i].t,bst[i].u,bst[i].v,bst[i].g);
  fclose(f);rename(tmp,ckptf);
}
static int load_ckpt(void){
  FILE*f=fopen(ckptf,"r");if(!f)return 0;
  ll bd,ev;int mm,bb,ff,nc;
  if(fscanf(f,"%lld %lld %d %d %d %d",&bd,&ev,&mm,&bb,&ff,&nc)!=6||mm!=m||bb!=B||ff!=F||nc>MAXC){fclose(f);return 0;}
  for(int i=0;i<nc;i++)if(fscanf(f,"%d %d %d %d",&bst[i].t,&bst[i].u,&bst[i].v,&bst[i].g)!=4){fclose(f);return 0;}
  fclose(f);
  nbst=nc;evals=ev;
  bstdef=evalcfg(bst,nbst);
  return 1;
}

static ll perturb_from(Cls*src,int nsrc,int extra){
  for(int att=0;att<60;att++){
    memcpy(cur,src,sizeof(Cls)*nsrc);ncur=nsrc;
    int k=2+rndi(3)+extra;
    for(int i2=0;i2<k;i2++){
      int ci=rndi(ncur);Cls*c=&cur[ci];
      if(c->t==LINK)c->g=rndi(m);
      else if(c->t==LOOP&&nlp)c->g=loopl[rndi(nlp)];
      else if(c->t==SEMI&&ninv)c->g=invl[rndi(ninv)];
    }
    if(rndi(100)<40)rand_struct();
    ll d=evalcfg(cur,ncur);
    if(d<INF)return d;
  }
  memcpy(cur,src,sizeof(Cls)*nsrc);ncur=nsrc;
  return evalcfg(cur,ncur);
}

static int getlist(Cls*c,const int**lst,int*buf){
  if(c->t==FIXO)return 0;
  if(c->t==LINK){for(int g=0;g<m;g++)buf[g]=g;*lst=buf;return m;}
  if(c->t==LOOP){*lst=loopl;return nlp;}
  *lst=invl;return ninv;
}
static ll dive3(ll d,double tend2){
  static int bi[MAXM],bj[MAXM],bk[MAXM];
  int imp=1;
  while(imp&&now()<tend2){
    imp=0;
    for(int i=0;i<ncur;i++){
      const int*li;int nli=getlist(&cur[i],&li,bi);
      if(!nli)continue;
      for(int j=i+1;j<ncur;j++){
        const int*lj;int nlj=getlist(&cur[j],&lj,bj);
        if(!nlj)continue;
        for(int k=j+1;k<ncur;k++){
          const int*lk;int nlk=getlist(&cur[k],&lk,bk);
          if(!nlk)continue;
          if((ll)nli*nlj*nlk>6000)continue;
          if(now()>tend2)return d;
          int ogi=cur[i].g,ogj=cur[j].g,ogk=cur[k].g;
          int bgi=ogi,bgj=ogj,bgk=ogk;ll bd=d;
          for(int a=0;a<nli;a++){cur[i].g=li[a];
            for(int b=0;b<nlj;b++){cur[j].g=lj[b];
              for(int c2=0;c2<nlk;c2++){
                if(li[a]==ogi&&lj[b]==ogj&&lk[c2]==ogk)continue;
                cur[k].g=lk[c2];
                ll dd=evalcfg(cur,ncur);
                if(dd<bd){bd=dd;bgi=li[a];bgj=lj[b];bgk=lk[c2];}
              }}}
          cur[i].g=bgi;cur[j].g=bgj;cur[k].g=bgk;
          if(bd<d){d=bd;imp=1;
            printf("dive3 improved to %lld evals=%lld\n",d,evals);fflush(stdout);
            if(d<bstdef){bstdef=d;memcpy(bst,cur,sizeof(Cls)*ncur);nbst=ncur;save_ckpt();}
          }
        }
      }
    }
  }
  return d;
}
static ll plateau_from(Cls*src,int nsrc){
  memcpy(cur,src,sizeof(Cls)*nsrc);ncur=nsrc;
  ll d=evalcfg(cur,ncur);
  int steps=6+rndi(12);
  for(int s=0;s<steps;s++){
    if(rndi(100)<30){
      Cls bak[MAXC];int nbak=ncur;memcpy(bak,cur,sizeof(Cls)*ncur);
      rand_struct();
      ll nd=evalcfg(cur,ncur);
      if(nd<=d)d=nd;
      else{memcpy(cur,bak,sizeof(Cls)*nbak);ncur=nbak;}
    }else{
      allow_side=1;
      ll nd=scan_class(cur,ncur,rndi(ncur),d);
      allow_side=0;
      if(nd<d)d=nd;
    }
  }
  return evalcfg(cur,ncur);
}
static int fixo_target(int j){
  for(int i=0;i<ncur;i++)if(cur[i].t==FIXO&&cur[i].u==j)return cur[i].v;
  return rndi(B);
}

static void dump_edges(void){
  ll d=evalcfg(bst,nbst);
  char ef[320];snprintf(ef,sizeof ef,"%s.edges",ckptf);
  FILE*f=fopen(ef,"w");if(!f)return;
  if(d<INF){
    for(int a=0;a<N;a++){
      u64 t=adjb[a][0];
      while(t){int b=__builtin_ctzll(t);t&=t-1;if(b>a)fprintf(f,"%d %d\n",a,b);}
      t=adjb[a][1];
      while(t){int b=64+__builtin_ctzll(t);t&=t-1;if(b>a)fprintf(f,"%d %d\n",a,b);}
    }
  }
  fclose(f);
  printf("DUMP def=%lld n=%d file=%s\n",d,N,ef);
}

int main(int argc,char**argv){
  if(argc<7){fprintf(stderr,"usage: lift GROUP B F SEED SECONDS CKPT\n");return 1;}
  rs=strtoull(argv[4],0,10)*2862933555777941757ULL+3037000493ULL;if(!rs)rs=1;
  snprintf(gname,sizeof gname,"%s",argv[1]);
  build_group(gname);
  B=atoi(argv[2]);F=atoi(argv[3]);
  double tlim=atof(argv[5]);
  snprintf(ckptf,sizeof ckptf,"%s",argv[6]);
  N=B*m+F;
  if(N>MAXN||B>MAXB){fprintf(stderr,"too big\n");return 1;}
  if(F>0&&m!=4){fprintf(stderr,"F>0 needs |G|=4\n");return 1;}
  nmask0=N>=64?~0ULL:((1ULL<<N)-1);
  nmask1=N<=64?0ULL:(N>=128?~0ULL:((1ULL<<(N-64))-1));
  int loaded=load_ckpt();
  if(tlim==0){
    if(!loaded){fprintf(stderr,"no ckpt to dump\n");return 1;}
    dump_edges();
    return 0;
  }
  if(tlim<0){
    if(!loaded){fprintf(stderr,"no ckpt to dive\n");return 1;}
    memcpy(cur,bst,sizeof(Cls)*nbst);ncur=nbst;
    ll d=evalcfg(cur,ncur);
    printf("dive start def=%lld\n",d);fflush(stdout);
    d=dive3(d,now()-tlim);
    if(d<bstdef){bstdef=d;memcpy(bst,cur,sizeof(Cls)*ncur);nbst=ncur;}
    save_ckpt();
    printf("RESULT group=%s B=%d F=%d N=%d best=%lld evals=%lld (dive)\n",gname,B,F,N,bstdef,evals);
    if(bstdef==0)dump_edges();
    return 0;
  }
  ll curdef;
  if(loaded){
    memcpy(cur,bst,sizeof(Cls)*nbst);ncur=nbst;
    curdef=evalcfg(cur,ncur);
    printf("resumed best=%lld evals=%lld nc=%d\n",bstdef,evals,ncur);
  }else{
    curdef=INF;
    for(int t2=0;t2<20000&&curdef>=INF;t2++){
      if(!init_struct())continue;
      curdef=evalcfg(cur,ncur);
    }
    if(curdef>=INF){fprintf(stderr,"init failed\n");return 1;}
    memcpy(bst,cur,sizeof(Cls)*ncur);nbst=ncur;bstdef=curdef;
    printf("init def=%lld nc=%d ninv=%d nlp=%d\n",curdef,ncur,ninv,nlp);
  }
  double t0=now(),tend=t0+tlim,lastck=t0;
  int stall=0,hops=0;
  while(bstdef>0){
    if(now()>=tend)break;
    ll tot=0;int nr=B+F;
    for(int ri=0;ri<nr;ri++)tot+=repdef[ri];
    int u;
    if(tot<=0||rndi(100)<20)u=rndi(B);
    else{
      ll pick=(ll)(rnd()%(u64)tot);int ri=0;
      while(pick>=repdef[ri]){pick-=repdef[ri];ri++;}
      u=ri<B?ri:fixo_target(ri-B);
    }
    int cand[MAXC],ncand=0;
    for(int i=0;i<ncur;i++){
      Cls*c=&cur[i];
      if(c->t==LINK&&(c->u==u||c->v==u))cand[ncand++]=i;
      else if((c->t==LOOP||c->t==SEMI)&&c->u==u)cand[ncand++]=i;
    }
    if(!ncand){stall++;continue;}
    int ci=cand[rndi(ncand)];
    allow_side=1;
    ll nd=scan_class(cur,ncur,ci,curdef);
    allow_side=0;
    if(nd<curdef){
      curdef=nd;evalcfg(cur,ncur);
      if(curdef<bstdef){
        bstdef=curdef;memcpy(bst,cur,sizeof(Cls)*ncur);nbst=ncur;
        printf("t=%.1f best=%lld evals=%lld\n",now()-t0,bstdef,evals);fflush(stdout);
        save_ckpt();hops=0;
      }
      stall=0;
    }else stall++;
    if(stall>30){
      curdef=polish(curdef,1);
      evalcfg(cur,ncur);
      if(curdef<bstdef){
        bstdef=curdef;memcpy(bst,cur,sizeof(Cls)*ncur);nbst=ncur;
        printf("t=%.1f best=%lld evals=%lld (polish)\n",now()-t0,bstdef,evals);fflush(stdout);
        save_ckpt();stall=0;hops=0;
      }else{
        hops++;
        if(hops>150){
          ll nd2=INF;
          for(int t2=0;t2<20000&&nd2>=INF;t2++){
            if(!init_struct())continue;
            nd2=evalcfg(cur,ncur);
          }
          if(nd2<INF)curdef=nd2;
          else curdef=perturb_from(bst,nbst,0);
          hops=0;
        }else if(bstdef<INF&&rndi(2)){
          curdef=plateau_from(bst,nbst);
        }else{
          int extra=hops/40;if(extra>4)extra=4;
          if(bstdef<INF&&rndi(100)<75)curdef=perturb_from(bst,nbst,extra);
          else curdef=perturb_from(cur,ncur,extra);
        }
        stall=0;
      }
    }
    if(now()-lastck>10){save_ckpt();lastck=now();}
  }
  save_ckpt();
  evalcfg(bst,nbst);
  printf("RESULT group=%s B=%d F=%d N=%d best=%lld evals=%lld repdef=",gname,B,F,N,bstdef,evals);
  for(int ri=0;ri<B+F;ri++)printf("%d,",repdef[ri]);
  printf("\n");
  if(bstdef==0)dump_edges();
  return 0;
}
