/* Exhaust u1,...,u5 on the flexible second-order subcase
   (u0,u6,u7,u8,u9)=(4,17,0,0,0). */

#include <stdio.h>

#define P 29

static const int equations[9][21] = {
 {25,8,1,1,7,8, 12,24,24,21,16, 10,18,24,2,27,0,15,6,25,14},
 {2,25,14,14,11,25, 23,17,17,4,21, 24,20,17,28,1,0,7,26,2,22},
 {21,20,18,20,21,5, 12,23,9,12,19, 0,10,5,12,18,14,23,8,28,4},
 {12,26,9,18,2,14, 18,5,7,18,4, 27,26,9,9,12,17,19,23,0,2},
 {3,12,17,0,20,17, 19,21,13,18,8, 28,9,0,17,16,4,3,1,23,26},
 {28,26,10,23,16,14, 4,14,12,11,16, 15,18,25,8,10,4,4,20,16,7},
 {21,21,21,13,9,5, 18,19,19,25,22, 25,20,13,16,0,21,17,13,18,22},
 {21,16,2,2,14,16, 24,19,19,13,3, 20,7,19,4,25,0,1,12,21,28},
 {19,20,17,17,3,20, 1,2,2,9,11, 25,16,2,5,24,0,23,15,19,6}
};

static int evaluate(int row, const int *u) {
    int cursor=0, total=equations[row][cursor++];
    for(int i=0;i<5;++i) total += equations[row][cursor++]*u[i];
    for(int i=0;i<5;++i) total += equations[row][cursor++]*u[i]*u[i];
    for(int i=0;i<5;++i)
        for(int j=i+1;j<5;++j)
            total += equations[row][cursor++]*u[i]*u[j];
    total %= P;
    return total < 0 ? total + P : total;
}

int main(void) {
    long tested=0, solutions=0;
    int u[5];
    for(u[0]=0;u[0]<P;++u[0])
    for(u[1]=0;u[1]<P;++u[1])
    for(u[2]=0;u[2]<P;++u[2])
    for(u[3]=0;u[3]<P;++u[3])
    for(u[4]=0;u[4]<P;++u[4]) {
        ++tested;
        int valid=1;
        for(int row=0;row<9;++row)
            if(evaluate(row,u)){valid=0;break;}
        if(valid){
            ++solutions;
            printf("solution u=[4,%d,%d,%d,%d,%d,17,0,0,0]\n",
                   u[0],u[1],u[2],u[3],u[4]);
        }
    }
    printf("tested=%ld solutions=%ld\n",tested,solutions);
    return solutions ? 0 : 2;
}
