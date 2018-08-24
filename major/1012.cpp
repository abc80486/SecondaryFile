#include"LinearEquationSolve.h"
#include<cstdio>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<time.h>

const double PI=3.14159246;
using namespace std;
double G[N][N]={0},B[N][N]={0},e[N],f[N],P[N]={0.0},Q[N]={0.0},dut[N],a[N],b[N];
double pline[N][N],qline[N][N];
void show2(int n,double (*a)[N]){
    //ofstream fout("out.txt");
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            
            if(j==0){cout<<"  "<<i<<"  ";}
            else{
                cout.width(10);
                if(i==0){cout<<j<<"  ";}
                else{
                //printf("%lf,",a[i][j]);
                cout<<a[i][j]<<"   ";}
            }
        }
        cout<<endl<<endl;
    }
}
void show1(int n,double *a){
    
    for(int j=1;j<=n;j++){
            //printf("%lf,",a[j]);
            cout.width(6);
            cout<<a[j]<<"   ";
     }
    //fout<<endl<<endl;
    cout<<endl<<endl;
}

void acquirepqab(int n,int m,double *lp,double *lq,double *pq){
    int k=1; 
    for(int i=1;i<=n+1+m;i++){
        lp[i]=lq[i]=0.0;
       
        a[i]=b[i]=0.0;
    }
   
    //show2(6,G);show2(6,B);
    //show1(6,e);show1(6,f);
    for(int i=2;i<=n+1;i++){
        for(int j=1;j<=n+1+m;j++){
     lp[i]+=e[i]*(G[i][j]*e[j]-B[i][j]*f[j])+f[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
     lq[i]+=f[i]*(G[i][j]*e[j]-B[i][j]*f[j])-e[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
        }
    } 
    for(int i=2;i<=n+1;i++){
        lp[i]=P[i]-lp[i];lq[i]=Q[i]-lq[i];
        pq[k++]=lp[i];
        pq[k++]=lq[i];
    }
    for(int i=n+2;i<=n+m+1;i++){
        for(int j=1;j<=n+m+1;j++){
     lp[i]+=e[i]*(G[i][j]*e[j]-B[i][j]*f[j])-f[i]*(G[i][j]*f[j]+B[i][j]*e[j]);
        }
        int du;
        du=dut[i-n-1]*dut[i-n-1]-e[i]*e[i]-f[i]*f[i];
        lp[i]=P[i]-lp[i];
        pq[k++]=lp[i];
        pq[k++]=du;
    }
    for(int i=2;i<=n+m+1;i++){
        for(int j=1;j<=n+m+1;j++){
            if(j!=i){
                a[i]+=G[i][j]*e[j]-B[i][j]*f[j];
               b[i]+=G[i][j]*f[j]+B[i][j]*e[j];
            }
        }
        a[i]+=G[i][i]*e[i]-B[i][i]*f[i];
        b[i]+=G[i][i]*f[i]+B[i][i]*e[i];
    }
    
}

void acquireJ(int n,int m,double (*H)[N],double (*N1)[N],double (*J)[N],double (*L)[N]){
    for(int i=2;i<=n+1;i++){
        for(int j=2;j<=n+m+1;j++){
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
    for(int i=n+2;i<=n+m+1;i++){
        for(int j=2;j<=n+m+1;j++){
           if(i!=j){
               //J[i][j]=-G[i][j]*e[i]+B[i][j]*f[i];
               H[i][j]=-B[i][j]*e[i]+G[i][j]*f[i];
                 N1[i][j]=G[i][j]*e[i]+B[i][i]*f[i];
                J[i][j]=0.0;L[i][j]=0.0;}
            else{
                 //J[i][j]=-G[i][j]*e[i]+B[i][j]*f[i]+a[i];
                H[i][j]=-B[i][j]*e[i]+G[i][j]*f[i]+b[i];
                 N1[i][j]=G[i][j]*e[i]+B[i][i]*f[i]+a[i];
                 J[i][i]=2.0*e[i];L[i][i]=2.0*f[i];
            }

        }
    }
    
}
void acquirejmat(int n,int m,double (*jj)[N],double (*H)[N],double (*N1)[N],double (*J)[N],double (*L)[N]){
    int k=1,p=1;
    for(int i=2;i<=n+m+1;i++){
        for(int j=2;j<=n+m+1;j++){
            jj[p][k++]=H[i][j];
            jj[p][k++]=N1[i][j];
        }
        k=1;p++;
        for(int j=2;j<=n+m+1;j++){
            jj[p][k++]=J[i][j];
            jj[p][k++]=L[i][j];
        }
        k=1;p++;
    }


}
void show3(int s){
     printf("\n电压P：");
    show1(s,e);
    printf("电压Q: ");
    show1(s,f);
    printf("电压大小:   ");
    for(int i=1;i<=s;i++){
        cout.width(10);
        cout<<sqrt(e[i]*e[i]+f[i]*f[i])<<"  ";
    }
    printf("\n\n电压相位(度):");
    for(int i=1;i<=s;i++){
                cout.width(10);
        cout<<atan2(f[i],e[i])*180.0/PI<<"  ";
    }
   for(int i=1;i<=s;i++){
       P[1]+=G[1][i]*e[i]-B[1][i]*f[i];
       Q[1]+=-(B[1][i]*e[i])-(G[1][i]*f[i]);
   }
   P[1]=P[1]*e[1]-Q[1]*f[1];
   Q[1]=P[1]*f[1]+Q[1]*e[1];
   cout<<endl<<endl<<"平衡节点功率：" <<P[1]<<"＋j"<<Q[1]<<endl;
   ///*
   for(int i=1;i<=s;i++){
       for(int j=1;j<=s;j++){
           if(i!=j){
               pline[i][j]=B[i][0]*f[i]-e[i]*G[i][0]-G[i][j]*(e[i]-e[j])-B[i][j]*(f[j]-f[i]);
               qline[i][j]=e[i]*B[i][0]+f[i]*G[i][0]+(e[i]-e[j])*B[i][j]-G[i][j]*(f[j]-f[i]);
               pline[i][j]=pline[i][j]*e[i]-qline[i][j]*f[i];
               qline[i][j]=pline[i][j]*f[i]+qline[i][j]*e[i];
           }
       }
   }
   cout<<endl<<"功率分布："<<endl;
   cout<<" "<<"有功功率P:"<<endl;
    show2(s,pline);
    cout<<endl;
    cout<<" "<<"无功功率Q:"<<endl;
    show2(s,qline);
    cout<<endl;
}
ifstream fin("102.txt");
void inputdata(int *n,int *m,double(*g)[N],double(*b)[N],double *p0,double *q0){
     
    int i,j,tnum=0,sp,sd=0;
    double y,u;
    fin>>*n>>*m>>tnum;
    sp=*n+*m+1;
    double k[sp][sp];
    int kp[sp][2];
    //输入变压器信息
    for(int t1=1;t1<=tnum;t1++){
           double kt;
           fin>>kp[t1][0]>>kp[t1][1]>>kt;
           k[kp[t1][0]][kp[t1][1]]=kt;
       }
    //形成导纳矩阵；
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
    //输入对地电容节点；
    fin>>sd;
    while(sd--){
        int i=1;
        double c0=0.0;
        fin>>i>>c0;
        b[i][0]-=c0;
    }
    //形成导纳矩阵对角线；
    for(int i=1;i<=sp;i++){
        for(int j=0;j<=sp;j++){
            if(i!=j){
                g[i][i]-=g[i][j];
                b[i][i]-=b[i][j];
            }
        }
    }
    //输入节点功率；
    for(int i=2;i<=*n+1;i++){
        fin>>p0[i]>>q0[i];
     }
     for(int i=*n+2;i<=sp;i++){
         fin>>p0[i]>>dut[i-*n-1];
     }
     //输入节点电压；
    for(int i=1;i<=sp;i++){
        fin>>e[i]>>f[i];}
}

int main(){
    double lp[N],lq[N],le[N],lf[N];
    double h[N][N],n1[N][N],j1[N][N],l[N][N];
    double jj[N][N],pq[N];
    //double pline[N][N],qline[N][N];//线路功率；
    double out[N]={0.0};
    int n,m,s,k=0,c=1,kp=1,kkl=0;
    clock_t uyt,uio;
    double hjk;
    uyt=clock();
    double *er,*ew;
    inputdata(&n,&m,G,B,P,Q);
    s=n+m+1;//show2(6,B);
    //show2(6,G);
    do
    {
    kp=1;c=1;
    acquirepqab(n,m,lp,lq,pq);//show1(10,pq);
    acquireJ(n,m,h,n1,j1,l);
    acquirejmat(n,m,jj,h,n1,j1,l);//show2(10,jj);
    output(2*(n+m),jj,pq,out);//show1(10,out); 
    //cout<<endl<<endl<<"节点电压"<<k<<": "<<endl;
    for(int i=1;i<=s;i++){
        //cout<<e[i]<<" ";
    }
     //cout<<endl;
    for(int i=1;i<=s;i++){
       // cout<<"+j"<<f[i]<<" ";
    }
    //cout<<endl<<endl<<"电压偏移"<<k<<": "<<endl;
    for(int i=2;i<=s;i++){
        //cout<<le[i]<<" ";
    }
    cout<<endl;
    for(int i=2;i<=s;i++){
        //cout<<"+j"<<lf[i]<<" ";
    }
   
    //cout<<endl<<endl;
    for(int i=2;i<=n+m+1;i++){
        lf[i]=out[c++];le[i]=out[c++];}
    for(int i=2;i<=n+m+1;i++){
        e[i]-=le[i]; f[i]-=lf[i];}
        //show1(6,e);
    er=max_element(out+1,out+2*(n+m));
    ew=min_element(out+1,out+2*(n+m));
    k++;
      }
    while(*er>1.0e-5||*ew<-1.0e-5);
    //while(kkl--);
    show3(s);
    uio=clock();
    hjk=(double)(uio-uyt)/CLOCKS_PER_SEC;
    cout<<hjk<<"s"<<endl;
    return 0;
}