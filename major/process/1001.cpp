#include<cstdio>
#include<iostream>
const int N=20;
double a[20][20];
using namespace std;
void zeroline(int j,int i,int n){
   double k;
    k=-a[j][i]/a[i][i];
    for(int t=i;t<=n;t++){
        a[j][t]+=a[i][t]*k;
    }
}
bool zero(int m,int n){
    bool abjustzero(int,int);
    for(int i=m;i<=n;i++){
        if(a[i][i]==0){
          if(abjustzero(i,n)==false){
            return false;}
        }
        for(int j=i+1;j<=n;j++){
           if(a[j][i]!=0) {
               zeroline(j,i,n);
           }
        }
    }
    return true;
}

bool abjustzero(int set,int end){
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
bool matrixch(int n){
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
bool matrixnubmul(int n,double b){
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            a[i][j]*=b;
        }
    }
    return true;
}
bool matrixmul(double (*p)[20],double (*y)[20],int a,int b,int c){
    int k=1;
    double s=0,r[20][20];
    for(int i=1;i<=a;i++){
        k=1;s=0.0;
        while(k<=c){
            for(int j=1;j<=b;j++){
                s+=p[i][j]*y[j][k];
            }
            r[i][k]=s;
            s=0.0;
            k++;
        }
    }
    for(int i=1;i<=a;i++){
            for(int j=1;j<=c;j++){
                printf("%lf  ",r[i][j]);
            }
            cout<<endl;
        } 
    return true;
}
 bool matrixback(double (*s)[N],int n){
     double p[N][N];
     for(int i=1;i<=n;i++){
         for(int j=1;j<=n;j++){
             p[i][j]=s[i][j];
         }
     }
     
 }
int main(){
    bool matrixch(int n);
    bool abjustzero(int,int);
    int n;
    double u[20][20];
    cin>>n;

    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            scanf("%lf",&a[i][j]);
        }
    }
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            scanf("%lf",&u[i][j]);
        }
    }
     matrixmul(u,a,2,2,2);
    //matrixch(n);
    double p=1.0;
    for(int i=1;i<=n;i++){
        p*=a[i][i];
    }
   //if(!zero(1,n)) p=0;
  // else{
        for(int i=1;i<=n;i++){
            for(int j=1;j<=n;j++){
                //printf("%lf  ",a[i][j]);
            }
            cout<<endl;
        } 
  // }
    //cout<<p<<endl;
    
    return 0;
}