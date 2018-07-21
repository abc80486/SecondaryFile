#include<cstdio>
#include<iostream>
#include<fstream>
const int N=20;
#define double long double
//double a[20][20];
using namespace std;
bool inversematrix(double (*s)[N],int n);
    void zeroline(int j,int i,int n,double(*a)[N]);
double determinant(int m,int n,double(*a)[N]);
  bool abjustzero(int set,int end,double(*a)[N]);
bool matrixch(int n,double(*a)[N]);
bool matrixnubmul(int n,double b,double(*a)[N]);
bool matrixmul(double (*p)[N],double (*y)[N],int a,int b,double *r);
double pc(const int n,int h,int k,double(*a)[N]);
void ax(int n,double(*a)[N]);

void show(int n,double(*p)[N]){
     for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            cout<<p[i][j];
        }
        cout<<endl;
    }
 }
void zeroline(int j,int i,int n,double(*a)[N]){
   double k;
    k=-a[j][i]/a[i][i];
    for(int t=i;t<=n;t++){
        a[j][t]+=a[i][t]*k;
    }
}
double determinant(int m,int n,double(*a)[N]){

    double p=1.0;
    if(((n-m)==1)){
        p=a[m][m]*a[n][n]-a[m][n]*a[n][m];
    }
    else{
    for(int i=m;i<=n;i++){
        if(a[i][i]==0){
          if(abjustzero(i,n,a)==false){
            return 0;}
        }
        for(int j=i+1;j<=n;j++){
           if(a[j][i]!=0) {
               zeroline(j,i,n,a);
           }
        }
    }
    for(int i=1;i<=n;i++){
        p*=a[i][i];
    }}
    return p;
}

bool abjustzero(int set,int end,double(*a)[N]){
    int q=0;
    if(a[set][set]==0){
         for(int t=1;t<=end;t++){
                if(a[set][t]!=0){
                    int p;
                    q=1;
                    for(int u=1;u<=end;u++){
                        p=a[u][set];
                        a[u][set]=-a[u][t];
                        a[u][t]=p;
                    }
                    break;
                }
            }
            if(q==0) return false;
    }
    return true;
}
bool matrixch(int n,double(*a)[N]){
    //int a[n];
    double p;
    for(int i=1;i<=n;i++){
        for(int j=i+1;j<=n;j++){
            p=a[i][j];
            a[i][j]=a[j][i];
            a[j][i]=p;
        }
    }
    return true;
}
bool matrixnubmul(int n,double b,double(*a)[N]){
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            a[i][j]*=b;
        }
    }
    return true;
}
bool matrixmul(double (*p)[N],double *y,int a,int b,double *r){
    double s=0.0;
    for(int i=1;i<=a;i++){
            s=0.0;
            for(int j=1;j<=b;j++){
                s+=p[i][j]*y[j];
            }
            r[i]=s;
    }
    for(int i=1;i<=a;i++){
       // printf("%lf  ",r[i]);
        //cout<<endl;
    } 
    return true;
}
 bool inversematrix(double (*s)[N],int n){
     double p[N][N],k;
     for(int i=1;i<=n;i++){
         for(int j=1;j<=n;j++){
             p[i][j]=s[i][j];
         }
     }
     k=determinant(1,n,p);
     //show(n,p);
     //cout<<k<<endl;
     if(k!=0.0){
         //matrixch(n,s);
        // show(n,s);
        ax(n,s);
        //show(n,s);
         matrixnubmul(n,1.0/k,s);
         //show(n,s);
     }
     else 
     {return false;}
     for(int i=1;i<=n;i++){
         for(int j=1;j<=n;j++){
            // p[i][j]=p[i][j];
             //cout<<s[i][j];
         }
         //cout<<endl;
     }
     return true;
 }
 

bool output(int n,double(*a)[N],double *u,double *out){
    if(!inversematrix(a,n)) {
        //cout<<"error!";
        return false;
    }
    else {
        matrixmul(a,u,n,n,out);
    }
    return true;
}
double pc(const int n,int h,int k,double(*a)[N]){
    double s[N][N];
    int p=0,b=0;
    for(int i=1;i<=n;i++){
        if(i!=h){ 
            p++;
            b=0;
            for(int j=1;j<=n;j++){
                //b=1;
                if(j!=k){
                    b=b+1;
                    s[p][b]=a[i][j];
                   
                }
            }
           
        }
    }
    //show(n-1,s);
    if(((h+k)%2)==0){
        return determinant(1,n-1,s);
        
    }
    else{
        return -determinant(1,n-1,s);
    }
}
void ax(int n,double(*a)[N]){
    double t[N][N];
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            t[j][i]=pc(n,i,j,a);
        }
    }
    for(int i=1;i<=n;i++){
         for(int j=1;j<=n;j++){
             a[i][j]=t[i][j];
         }
     }
}

int main(){
    ifstream fin("jz.txt");
    int n;
    double a[N][N];
    double u[N];
    double out[N];
    
    fin>>n;
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            fin>>a[i][j];
        }
    }
    //cout<<endl;
    for(int i=1;i<=n;i++){
         fin>>u[i];
        
    }
    output(n,a,u,out);
    for(int i=1;i<=n;i++){
        cout<<out[i]<<endl;
    }
    //show(n,out);
    return 0;
}