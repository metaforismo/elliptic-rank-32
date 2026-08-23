/* Exhaust the reduced second-order obstruction system for P3 #1.

   The certified Groebner basis forces u6=21 and u5=25.  After substitution,
   six nonzero quadrics remain in u0,...,u4.  The fourth is affine-linear in
   u0, so only 29^4 outer tuples and a small exceptional branch are required.
*/

#include <stdio.h>

#define P 29

static const int equations[6][21] = {
    {27,10,13,0,0,4, 0,17,19,23,15, 15,9,16,18,12,15,8,15,13,6},
    {6,8,16,23,27,12, 0,14,2,26,13, 9,16,19,22,25,26,3,13,15,7},
    {24,14,25,12,23,27, 0,18,12,24,6, 20,16,24,9,16,11,20,14,0,27},
    {19,22,28,13,6,2, 0,14,14,7,11, 0,3,5,2,28,24,21,4,27,8},
    {24,0,3,20,11,4, 0,18,5,2,10, 8,3,18,11,15,21,6,3,11,25},
    {4,2,27,0,15,20, 1,4,3,8,2, 4,20,8,25,10,10,16,27,5,11}
};

static int mod(int value) {
    value %= P;
    return value < 0 ? value + P : value;
}

static void monomials(int u0, int u1, int u2, int u3, int u4, int *m) {
    int u[5] = {u0,u1,u2,u3,u4};
    int cursor = 0;
    m[cursor++] = 1;
    for (int i=0;i<5;++i) m[cursor++] = u[i];
    for (int i=0;i<5;++i) m[cursor++] = u[i]*u[i] % P;
    for (int i=0;i<5;++i)
        for (int j=i+1;j<5;++j)
            m[cursor++] = u[i]*u[j] % P;
}

static int evaluate(int row, const int *m) {
    int total = 0;
    for (int i=0;i<21;++i) total += equations[row][i]*m[i];
    return mod(total);
}

int main(void) {
    int inverse[P] = {0};
    for (int value=1; value<P; ++value)
        for (int candidate=1; candidate<P; ++candidate)
            if (value*candidate % P == 1) inverse[value] = candidate;
    long outer = 0, exceptional_u0 = 0, tested = 0, solutions = 0;
    for (int u1=0;u1<P;++u1)
    for (int u2=0;u2<P;++u2)
    for (int u3=0;u3<P;++u3)
    for (int u4=0;u4<P;++u4) {
        ++outer;
        /* Equation row 3: denominator is the coefficient of u0. */
        int denominator = mod(22 + 3*u2 + 5*u3 + 2*u4);
        int m0[21];
        monomials(0,u1,u2,u3,u4,m0);
        int constant = evaluate(3,m0);
        int first = 0, last = P;
        if (denominator) {
            first = mod(-constant * inverse[denominator]);
            last = first + 1;
        } else if (constant) {
            continue;
        } else {
            exceptional_u0 += P;
        }
        for (int u0=first; u0<last; ++u0) {
            int m[21];
            ++tested;
            monomials(u0,u1,u2,u3,u4,m);
            int valid = 1;
            for (int row=0;row<6;++row)
                if (evaluate(row,m)) { valid=0; break; }
            if (valid) {
                ++solutions;
                printf("solution u=[%d,%d,%d,%d,%d,25,21]\n",u0,u1,u2,u3,u4);
            }
        }
    }
    printf("outer=%ld exceptional_u0_cases=%ld tested=%ld solutions=%ld\n",
           outer, exceptional_u0, tested, solutions);
    return solutions ? 0 : 2;
}
