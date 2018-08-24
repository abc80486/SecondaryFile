#include<iostream>
#include<cstdio>
//#include<cstring>
#include<fstream>
//#include<math.h>
//#define double float
const int N=20;
//#define M 100
using namespace std;
double G[N][N]={0},B[N][N]={0},e[N],f[N],P[N],Q[N],a[N],b[N];
void show2(int n,double (*a)[N]){
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            printf("%lf ",a[i][j]);
        }
        cout<<endl<<endl;
    }
}
void show1(int n,double *a){
    for(int j=1;j<=n;j++){
            printf("%lf ",a[j]);
     }
    cout<<endl<<endl;
}

void inputdata(int *n,int *m,double(*g)[N],double(*b)[N],double *p0,double *q0){
    ifstream fin("100a.txt");
    int i,j,tnum=0,sp;
    double y,u;
    fin>>*n>>*m>>tnum;
    sp=*n+*m+1;
    double k[sp][sp];
    int kp[sp][2];
    for(int t1=1;t1<=tnum;t1++){
           double kt;
           fin>>kp[t1][0]>>kp[t1][1]>>kt;
           k[kp[t1][0]][kp[t1][1]]=kt;
       }
    while(fin>>i>>j&& i!=0 && j!=0){
       fin>>y>>u;
        g[i][j]=-y/(y*y+u*u);
        b[i][j]=u/(y*y+u*u);
       for(int ko=1;ko<=tnum;ko++){
           if(i==kp[ko][0]&&j==kp[ko][1]){
               g[i][0]=g[i][j]*(1.0-k[i][j])/k[i][j]/k[i][j];
               b[i][0]=b[i][j]*(1.0-k[i][j])/k[i][j]/k[i][j];
               g[j][0]=g[i][j]*(k[i][j]-1.0)/k[i][j];
               b[j][0]=b[i][j]*(k[i][j]-1.0)/k[i][j];
               g[i][j]/=k[i][j];
               b[i][j]/=k[i][j];
               break;
           }
           if(j==kp[ko][0]&&i==kp[ko][1]){
               g[j][0]=g[i][j]*(1-k[j][i])/k[j][i]/k[j][i];
               b[j][0]=b[i][j]*(1-k[j][i])/k[j][i]/k[j][i];
               g[i][0]=g[i][j]*(k[j][i]-1)/k[j][i];
               b[i][0]=b[i][j]*(k[j][i]-1)/k[j][i];
               g[i][j]/=k[j][i];
               b[i][j]/=k[j][i];
            break;
           }
       }
       g[j][i]=g[i][j];
       b[j][i]=b[i][j];
    }
    for(int i=1;i<=sp;i++){
        for(int j=0;j<=sp;j++){
            if(i!=j){
                
                g[i][i]-=g[i][j];
                b[i][i]-=b[i][j];
                
            }
        }
    }
    for(int i=2;i<=*n+1;i++){
        fin>>p0[i]>>q0[i];}
    for(int i=1;i<=*m;i++){
        fin>>p0[*n+1+i]>>q0[*n+1+i];
    }
}
void input( int n, int m,double *e,double *f){
    ifstream fin("100b.txt");
    int k=m+n+1;
    for(int i=1;i<=n+m+1;i++){
        fin>>e[i]>>f[i];}
    
    
}
void acquirepqab(int n,int m,double *lp,double *lq){
    for(int i=1;i<=n+1;i++){
        lp[i]=lq[i]=0.0;
    }
    for(int i=2;i<=n+1;i++){
        for(int j=1;j<=n+1;j++){
     lp[i]+=e[i]*(G[i][j]*e[j]-B[i][j]*f[j])+f[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
     lq[i]+=f[i]*(G[i][j]*e[j]-B[i][j]*f[j])-e[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
        }
    }
    for(int i=2;i<=n+1;i++){
        a[i]=(lp[i]*e[i]+(-lq[i])*(-1.0*f[i]))/(e[i]*e[i]+f[i]*f[i]);
        b[i]=((-1.0*lq[i])*e[i]-lp[i]*(-1.0*f[i]))/(e[i]*e[i]+f[i]*f[i]);
       lp[i]=P[i]-lp[i];lq[i]=Q[i]-lq[i];
    }
}

void acquireJ(int n,double (*H)[N],double (*N1)[N],double (*J)[N],double (*L)[N]){
    for(int i=2;i<=n+1;i++){
        for(int j=2;j<=n+1;j++){
            if(i==j){
                H[i][j]=-B[i][j]*e[i]+G[i][j]*f[i]+b[i];
                 N1[i][j]=G[i][j]*e[i]+B[i][i]*f[i]+a[i];
                 J[i][j]=-G[i][j]*e[i]+B[i][j]*f[i]+a[i];
                L[i][j]=-B[i][j]*e[i]+G[i][j]*f[i]-b[i];
            }
            else{
                H[i][j]=-B[i][j]*e[i]+G[i][j]*f[i];
                 N1[i][j]=G[i][j]*e[i]+B[i][i]*f[i];
                 J[i][j]=-G[i][j]*e[i]+B[i][j]*f[i];
                L[i][j]=-B[i][j]*e[i]+G[i][j]*f[i];
            }
     }
    }
}
void acquirejmat(int n,double (*jj)[N],double (*H)[N],double (*N1)[N],double (*J)[N],double (*L)[N]){
    int k=1,p=1;
    for(int i=2;i<=n+1;i++){
        for(int j=2;j<=n+1;j++){
            jj[p][k++]=H[i][j];
            jj[p][k++]=N1[i][j];
        }
        k=1;p++;
        for(int j=2;j<=n+1;j++){
            jj[p][k++]=J[i][j];
            jj[p][k++]=L[i][j];
        }
        k=1;p++;
    }


}
int main(){
    //double lp[N],lq[N];
    //double h[N][N],n1[N][N],j[N][N],l[N][N];
   
    bool createj(double (*jj)[N],double *pqu,int *n,int *m);
    double jj[N][N];
    double pq[N];
    int n,m;
    createj(jj,pq,&n,&m);
    //show1(n*2,pq); 
    //ofstream fout("100b.txt");
    
    //s=n+m+1;
    //show2(s,G);show2(s,B);show1(s,e);show1(s,f);
    //show1(2*n,pq);
    return 0;
}
bool createj(double (*jj)[N],double *pqu,int *n,int *m){
    double lp[N],lq[N];
    double h[N][N],n1[N][N],j[N][N],l[N][N];
    int s,k=1;
    int nn,mm;
    inputdata(&nn,&mm,G,B,P,Q);
    *n=nn;*m=mm;
    s=*m+*n+1;
    input(nn,mm,e,f); 
    acquirepqab(nn,mm,lp,lq);
    acquireJ(nn,h,n1,j,l);
    acquirejmat(nn,jj,h,n1,j,l);//show2(2*nn,jj);
    for(int i=2;i<=nn+1;i++){
            pqu[k++]=lp[i];
            pqu[k++]=lq[i];
     }
     //show1(2*nn,pqu);
     return true;
}

