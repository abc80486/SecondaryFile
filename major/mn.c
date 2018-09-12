
#include<stdio.h>
int main(){
    int a[3][2];
    int i,max=1,s=0;
    for(i=0;i<3;i++){
        scanf("%d%d",&a[i][0],&a[i][1]);
        max*=a[i][0];
    }
    int k;
    for(i=0;i<3;i++){
        k=max/a[i][0];
        while((k%a[i][0])!=a[i][1]) k+=(max/a[i][0]);
        s+=k;
    }
    printf("%d\n",s%max);
    return 0;
}
/*
注意：输入中数组第一列“3 5 7”中的任意两个元素的最大公约数必须为1；
第二列不能出现0并且模数合法；
如合法例子：       错误例子：因为4和6的最大公约数为2；
3 2    7 3        4 2
5 3    8 5       6 3
7 2    9 6       7 5  
23     213
本题还可以推广，即模数不一定为3，如：
4          //代表组数；
3 2
5 3
7 2
8 6  //答案758；
*/
/*
#include<stdio.h>
int main(){
    int n,i,max=1,s=0;
    int a[1000][2];
    scanf("%d",&n);
    for(i=0;i<n;i++){
        scanf("%d%d",&a[i][0],&a[i][1]);
        max*=a[i][0];
    }
    int k;
    for(i=0;i<n;i++){
        k=max/a[i][0];
        while((k%a[i][0])!=a[i][1]) k+=(max/a[i][0]);
        s+=k;
    }
    printf("%d\n",s%max);
    return 0;
}
*/



