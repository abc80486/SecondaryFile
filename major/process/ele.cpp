#include<iostream>
#include<cstdio>
#include<cstring>
#include<fstream>
#include<math.h>

#define N 100
#define M 100
using namespace std;
 double G[N][N],B[N][N],e[N],f[N],P[M],Q[M],u1=1.06;;
int note,pqnote;
void inputdata(int *n,int *m,double(*g)[N],double(*b)[N]){
    int i,j;
    double y,u;
    scanf("%d%d",n,m);
    while(scanf("%d,%d",&i,&j)==2 && i!=0 && j!=0){
        scanf("%lf %lf",&y,&u);
        g[i][j]=y/sqrt(y*y+u*u);
        b[i][j]=-u/sqrt(y*y+u*u);
        if(i>j) *n=i;
        else *n=j;
    }
}
void input( int n, int m,double *e,double *f,double *p0,double *q0){
    int k=m+n+1;
    while(--k){
        scanf("%lf%lf",&e[k],&f[k]);}
    while(--n){
        scanf("%lf%lf",&p0[n+1],&q0[n+1]);}
    
}
void acquirepq(int n,int m,double *lp,double *lq){
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n+m+1;j++){
     lp[i]+=e[i]*(G[i][j]*e[j]-B[i][j]*f[j])+f[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
     lq[i]+=f[i]*(G[i][j]*e[j]-B[i][j]*f[j])-e[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
        }
    }
    for(int i=1;i<=n;i++){
       lp[i]=P[i]-lp[i];lq[i]=Q[i]-lq[i];
    }
}
void acquireab(int n,double *a,double *b){
    for(int i=1;i<=n;i++){
        a[i]=(P[i]*e[i]+(-Q[i])*(-f[i]))/(e[i]*e[i]+f[i]*f[i]);
        b[i]=((-Q[i])*e[i]-P[i]*(-f[i]))/(e[i]*e[i]+f[i]*f[i]);
    }

}

void DataInput(){
    ifstream fin("indata.txt");
    fin>>note>>pqnote;
    for(int j=1;j<=note;j++)
        for(int i=1;i<=note;i++) fin>>G[j][i];
    for(int j=1;j<=note;j++)
        for(int i=1;i<=note;i++) fin>>B[j][i];
    for(int i=1;i<=note;i++) fin>>e[i];
        
    for(int i=1;i<=note;i++) fin>>f[i];

    for(int j=1;j<=pqnote;j++) fin>>P[j]>>Q[j];
 
}
int main(){
     double  P0[M]={0.0},Q0[M]={0.0},DP[M],DQ[M],aii[M],bii[M];
    DataInput();

    for(int i=2;i<=note;i++){
        for(int j=1;j<=note;j++)    {
    P0[i]+=e[i]*(G[i][j]*e[j]-B[i][j]*f[j])+f[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
     Q0[i]+=f[i]*(G[i][j]*e[j]-B[i][j]*f[j])-e[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
        }
    }
    for(int i=2;i<=note;i++){
        DP[i]=P[i]-P0[i];
        DQ[i]=Q[i]-Q0[i];
    }
    for(int i=2;i<=note;i++){
        aii[i]=(P0[i]*e[i]+(-Q0[i])*(-f[i]))/(e[i]*e[i]+f[i]*f[i]);
        bii[i]=((-Q0[i])*e[i]-P0[i]*(-f[i]))/(e[i]*e[i]+f[i]*f[i]);
    }
    for(int i=2;i<=note;i++) printf("%f ",P0[i]); 
    printf("\n"); 
    while(1);
    return 0;
   
}