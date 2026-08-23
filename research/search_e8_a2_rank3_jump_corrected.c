/* Exact finite-field search for the corrected E8+A2^3 section ansatz.

       x = q^2 + r,
       y = q^3 + (3/2) q r + h,

   where deg q=2 and deg r,deg h<=1.  The (3/2)qr term is essential: without
   it the coefficients of t^9 and t^8 force r=0 on q2!=0.

   At t=0 put x(0)=u^2,y(0)=u^3.  This gives

       r0=u^2-q0^2,
       h0=(u-q0)^2(2u+q0)/2.

   For u!=0 the t coefficient solves

       r1=q1(u-q0)+1/(3 q2^3 (u-q0));

   u=q0 is impossible.  For u=0 the t coefficient vanishes identically and
   r1 remains free.  The t^7 and t^6 coefficients solve h1 and mu.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define WIDTH 13

static int P;
static int mod_i64(int64_t x) { int r=(int)(x%P); return r<0?r+P:r; }
static int addm(int a,int b){return mod_i64((int64_t)a+b);}
static int subm(int a,int b){return mod_i64((int64_t)a-b);}
static int mulm(int a,int b){return mod_i64((int64_t)a*b);}
static int powm(int a,int n){int r=1;while(n){if(n&1)r=mulm(r,a);a=mulm(a,a);n>>=1;}return r;}
static int invm(int a){return a?powm(a,P-2):0;}
static void zero(int *a){for(int i=0;i<WIDTH;i++)a[i]=0;}
static void pmul(const int *a,const int *b,int *c){
    zero(c);
    for(int i=0;i<WIDTH;i++) if(a[i])
      for(int j=0;i+j<WIDTH;j++) if(b[j])
        c[i+j]=addm(c[i+j],mulm(a[i],b[j]));
}

static int verify(
    int q0,int q1,int q2,int r0,int r1,int h0,int h1,int lam,int mu
){
    int q[WIDTH],r[WIDTH],h[WIDTH],q2p[WIDTH],q3p[WIDTH],qr[WIDTH];
    int x[WIDTH],y[WIDTH],x2p[WIDTH],x3p[WIDTH],y2p[WIDTH];
    int a[WIDTH],b[WIDTH],c[WIDTH],f[WIDTH],fsq[WIDTH];
    zero(q);zero(r);zero(h);
    q[0]=q0;q[1]=q1;q[2]=q2;r[0]=r0;r[1]=r1;h[0]=h0;h[1]=h1;
    pmul(q,q,q2p);pmul(q2p,q,q3p);pmul(q,r,qr);
    int three_halves=mulm(3,invm(2));
    for(int i=0;i<WIDTH;i++){
        x[i]=addm(q2p[i],r[i]);
        y[i]=addm(addm(q3p[i],mulm(three_halves,qr[i])),h[i]);
    }
    pmul(x,x,x2p);pmul(x2p,x,x3p);pmul(y,y,y2p);
    zero(a);zero(b);zero(c);a[1]=1;b[0]=P-1;b[1]=1;c[0]=subm(0,lam);c[1]=1;
    pmul(a,b,f);pmul(f,c,a);pmul(a,a,fsq);
    zero(b);b[0]=subm(0,mu);b[1]=1;pmul(fsq,b,f);
    for(int i=0;i<WIDTH;i++) if(subm(subm(y2p[i],x3p[i]),f[i])) return 0;
    return 1;
}

static int compute_h1(int q2){return mulm(invm(2),invm(mulm(mulm(q2,q2),q2)));}
static int compute_mu(int lam,int q1,int q2,int r1,int h0,int h1){
    int q22=mulm(q2,q2),q23=mulm(q22,q2),mu=0;
    mu=subm(mu,mulm(2,mulm(h0,q23)));
    mu=subm(mu,mulm(6,mulm(mulm(h1,q1),q22)));
    mu=subm(mu,mulm(2,lam));
    mu=addm(mu,mulm(mulm(3,invm(4)),mulm(q22,mulm(r1,r1))));
    mu=subm(mu,2);
    return mu;
}

int main(int argc,char **argv){
    if(argc<2||argc>4){fprintf(stderr,"usage: %s PRIME [MAX_PRINT] [smooth|cusp|both]\n",argv[0]);return 2;}
    P=atoi(argv[1]);int max_print=argc>=3?atoi(argv[2]):100;
    const char *mode=argc==4?argv[3]:"both";
    int smooth=mode[0]=='s'||mode[0]=='b',cusp=mode[0]=='c'||mode[0]=='b';
    if(P<=3||P==79||(!smooth&&!cusp)){fprintf(stderr,"bad prime or mode\n");return 2;}
    uint64_t sn=0,sg=0,cn=0,cg=0,so=0,co=0;
    int printed=0,inv2=invm(2),inv3=invm(3);
    if(smooth) for(int q0=0;q0<P;q0++)for(int q1=0;q1<P;q1++)for(int q2=1;q2<P;q2++){
      int q23=mulm(mulm(q2,q2),q2),h1=compute_h1(q2);
      for(int u=1;u<P;u++){
        so++;int d=subm(u,q0);if(!d)continue;
        int r0=subm(mulm(u,u),mulm(q0,q0));
        int h0=mulm(inv2,mulm(mulm(d,d),addm(mulm(2,u),q0)));
        int r1=addm(mulm(q1,d),mulm(inv3,invm(mulm(q23,d))));
        for(int lam=0;lam<P;lam++){
          int mu=compute_mu(lam,q1,q2,r1,h0,h1);
          if(!verify(q0,q1,q2,r0,r1,h0,h1,lam,mu))continue;
          sn++;int generic=lam!=0&&lam!=1&&mu!=0&&mu!=1&&mu!=lam;if(generic)sg++;
          if(printed<max_print){
            printf("R3CPOINT|branch=smooth|q0=%d|q1=%d|q2=%d|u=%d|r0=%d|r1=%d|h0=%d|h1=%d|lambda=%d|mu=%d|generic=%d\n",q0,q1,q2,u,r0,r1,h0,h1,lam,mu,generic);printed++;
          }
        }
      }
    }
    if(cusp) for(int q0=0;q0<P;q0++)for(int q1=0;q1<P;q1++)for(int q2=1;q2<P;q2++){
      int h1=compute_h1(q2),r0=subm(0,mulm(q0,q0));
      int h0=mulm(inv2,mulm(mulm(q0,q0),q0));
      for(int r1=0;r1<P;r1++){
        co++;
        for(int lam=0;lam<P;lam++){
          int mu=compute_mu(lam,q1,q2,r1,h0,h1);
          if(!verify(q0,q1,q2,r0,r1,h0,h1,lam,mu))continue;
          cn++;int generic=lam!=0&&lam!=1&&mu!=0&&mu!=1&&mu!=lam;if(generic)cg++;
          if(printed<max_print){
            printf("R3CPOINT|branch=cusp|q0=%d|q1=%d|q2=%d|u=0|r0=%d|r1=%d|h0=%d|h1=%d|lambda=%d|mu=%d|generic=%d\n",q0,q1,q2,r0,r1,h0,h1,lam,mu,generic);printed++;
          }
        }
      }
    }
    printf("R3CSUMMARY|p=%d|smooth_outer=%llu|smooth_identities=%llu|smooth_generic=%llu|cusp_outer=%llu|cusp_identities=%llu|cusp_generic=%llu|max_print=%d\n",P,(unsigned long long)so,(unsigned long long)sn,(unsigned long long)sg,(unsigned long long)co,(unsigned long long)cn,(unsigned long long)cg,max_print);
    return 0;
}
