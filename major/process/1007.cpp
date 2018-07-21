#include"LinearEquationSolve.h"
#include"createjmatrix.h"
#include<cstdio>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<cmath>
const double PI=3.14159246;
using namespace std;
int main(){
    ifstream fin("100b.txt");
    double lp[N],lq[N],le[N],lf[N];
    double jj[N][N],pq[N],e[N],f[N];
    double pline[N][N],qline[N][N];//线路功率；
    double out[N]={0.0};
    int n,m,s,k=0,c=1;;
    double *er,*ew;
    //show2(s,G);
    do
    {
        for(int i=1;i<=5;i++){
            fin>>e[i]>>f[i];
        }
    createj(jj,pq,&n,&m,e,f);
    //show2(2*n+2*m,jj);
    cout<<ends;
    output(2*n,jj,pq,out);
    show1(n,e);
    c=1;
    for(int i=2;i<=n+m;i++){
        lf[i]=out[c++];
        le[i]=out[c++];
     }
     //show1(n,)
     //ofstream fout("100b.txt");
    for(int i=2;i<=n+m;i++){
        e[i]+=le[i];
        f[i]+=lf[i];
        //fout<<e[i]<<" "<<f[i]<<endl;
    }
    er=max_element(out+1,out+2*n);
    ew=min_element(out+1,out+2*n);
    k++;
      }
    while(*er>1.0e-5||*ew<-1.0e-5);

   // printf("电压P：");
   // show1(s,e);
    //printf("电压Q: ");
    //show1(s,f);
    //printf("电压大小:   ");
    for(int i=1;i<=s;i++){
        //cout<<sqrt(e[i]*e[i]+f[i]*f[i])<<" ";
    }
    //printf("\n\n电压相位(度):");
    for(int i=1;i<=s;i++){
        //cout<<atan2(f[i],e[i])*180.0/PI<<" ";
    }
    cout<<endl;
   //show2(s,B);
  //cout<<endl<<P[s]<<" "<<Q[s];
   for(int i=1;i<=s;i++){
       P[s]+=G[s][i]*e[i]-B[s][i]*f[i];
       Q[s]+=-(B[s][i]*e[i])-(G[s][i]*f[i]);
       //cout<<B[s][i]<<" "<<G[s][i]<<endl;
      // cout<<Q[s]<<endl;
   }
   //Q[s]=Q[s]*e[s];
  // cout<<endl<<P[s]<<" "<<Q[s];
   P[s]=P[s]*e[s]-Q[s]*f[s];
   Q[s]=P[s]*f[s]+Q[s]*e[s];
   
  // cout<<endl<<endl<<"平衡节点功率：" <<P[s]<<"＋j"<<Q[s]<<endl;
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
    //show2(s,qline);
   //*/
    cout<<endl;
    return 0;
}