/* Exhaust smooth/cuspidal t=0 branches of the E8+A2^3 rank-jump system.

   The section ansatz is x=q^2+r, y=q^3+s.  Write its smooth specialization
   on the cuspidal fiber t=0 as x(0)=u^2, y(0)=u^3 with u != 0.  Then

       r0 = u^2-q0^2,  s0 = u^3-q0^3,

   and the t^1 residual is affine-linear in r1.  The t^7 and t^6 residuals
   solve for s1 and mu.  We enumerate the remaining q0,q1,q2,u,lambda exactly
   over GF(p) and verify the complete polynomial identity coefficientwise.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define WIDTH 13

static int P;

static int mod_i64(int64_t value) {
    int result = (int)(value % P);
    return result < 0 ? result + P : result;
}

static int add(int a, int b) { return mod_i64((int64_t)a + b); }
static int sub(int a, int b) { return mod_i64((int64_t)a - b); }
static int mul(int a, int b) { return mod_i64((int64_t)a * b); }

static int power(int base, int exponent) {
    int result = 1;
    while (exponent) {
        if (exponent & 1) result = mul(result, base);
        base = mul(base, base);
        exponent >>= 1;
    }
    return result;
}

static int inverse(int value) {
    if (!value) return 0;
    return power(value, P - 2);
}

static void clear_poly(int *poly) {
    for (int i = 0; i < WIDTH; ++i) poly[i] = 0;
}

static void poly_mul(const int *left, const int *right, int *result) {
    clear_poly(result);
    for (int i = 0; i < WIDTH; ++i) {
        if (!left[i]) continue;
        for (int j = 0; i + j < WIDTH; ++j) {
            if (!right[j]) continue;
            result[i + j] = add(result[i + j], mul(left[i], right[j]));
        }
    }
}

static int verify_identity(
    int q0, int q1, int q2, int r0, int r1, int s0, int s1, int lam, int mu
) {
    int q[WIDTH], r[WIDTH], s[WIDTH], q_sq[WIDTH], q_cube[WIDTH];
    int x[WIDTH], y[WIDTH], x_sq[WIDTH], x_cube[WIDTH], y_sq[WIDTH];
    int f0[WIDTH], f1[WIDTH], f2[WIDTH], fiber[WIDTH], fiber_sq[WIDTH];
    clear_poly(q); clear_poly(r); clear_poly(s);
    q[0] = q0; q[1] = q1; q[2] = q2;
    r[0] = r0; r[1] = r1;
    s[0] = s0; s[1] = s1;
    poly_mul(q, q, q_sq);
    poly_mul(q_sq, q, q_cube);
    for (int i = 0; i < WIDTH; ++i) {
        x[i] = add(q_sq[i], r[i]);
        y[i] = add(q_cube[i], s[i]);
    }
    poly_mul(x, x, x_sq);
    poly_mul(x_sq, x, x_cube);
    poly_mul(y, y, y_sq);

    clear_poly(f0); clear_poly(f1); clear_poly(f2);
    f0[0] = 0; f0[1] = 1;
    f1[0] = P - 1; f1[1] = 1;
    f2[0] = sub(0, lam); f2[1] = 1;
    poly_mul(f0, f1, fiber);
    poly_mul(fiber, f2, f0);
    poly_mul(f0, f0, fiber_sq);
    clear_poly(f1); f1[0] = sub(0, mu); f1[1] = 1;
    poly_mul(fiber_sq, f1, fiber);
    for (int i = 0; i < WIDTH; ++i) {
        if (sub(sub(y_sq[i], x_cube[i]), fiber[i]) != 0) return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "usage: %s PRIME [MAX_PRINT] [smooth|cusp|both]\n", argv[0]);
        return 2;
    }
    P = atoi(argv[1]);
    int max_print = argc >= 3 ? atoi(argv[2]) : 100;
    const char *branch = argc == 4 ? argv[3] : "both";
    int run_smooth = branch[0] == 's' || branch[0] == 'b';
    int run_cusp = branch[0] == 'c' || branch[0] == 'b';
    if (!run_smooth && !run_cusp) {
        fprintf(stderr, "branch must be smooth, cusp, or both\n");
        return 2;
    }
    if (P <= 3 || P == 79) {
        fprintf(stderr, "choose a good odd prime other than 3 and 79\n");
        return 2;
    }
    int inv2 = inverse(2);
    uint64_t outer = 0, denominator_zero = 0, compatible_zero = 0, identities = 0;
    uint64_t generic_identities = 0;
    if (run_smooth) for (int q0 = 0; q0 < P; ++q0) {
        for (int q1 = 0; q1 < P; ++q1) {
            for (int q2 = 1; q2 < P; ++q2) {
                int q2_2 = mul(q2, q2), q2_3 = mul(q2_2, q2);
                for (int u = 1; u < P; ++u) {
                    ++outer;
                    int r0 = sub(mul(u, u), mul(q0, q0));
                    int s0 = sub(mul(mul(u, u), u), mul(mul(q0, q0), q0));
                    int constant = 0;
                    constant = add(constant, mul(6, mul(mul(mul(q0, q0), q1), q2_3)));
                    constant = add(constant, mul(6, mul(mul(mul(q0, q1), q2_3), u)));
                    constant = sub(constant, mul(12, mul(mul(mul(q1, q2_3), u), u)));
                    constant = sub(constant, 1);
                    int coefficient = 0;
                    coefficient = sub(coefficient, mul(12, mul(q0, q2_3)));
                    coefficient = sub(coefficient, mul(18, mul(mul(q1, q1), q2_2)));
                    coefficient = add(coefficient, mul(3, mul(q2_3, u)));
                    int r1_start, r1_end;
                    if (!coefficient) {
                        ++denominator_zero;
                        if (constant) continue;
                        ++compatible_zero;
                        r1_start = 0;
                        r1_end = P;
                    } else {
                        r1_start = mul(sub(0, constant), inverse(coefficient));
                        r1_end = r1_start + 1;
                    }
                    for (int r1_raw = r1_start; r1_raw < r1_end; ++r1_raw) {
                      int r1 = r1_raw % P;
                      int s1_numerator = 0;
                      s1_numerator = add(s1_numerator, mul(12, mul(mul(mul(q1, q2_3), r0), 1)));
                      s1_numerator = add(s1_numerator, mul(18, mul(mul(mul(mul(q1, q1), q2_2), r1), 1)));
                      s1_numerator = add(s1_numerator, mul(12, mul(mul(mul(q0, q2_3), r1), 1)));
                      s1_numerator = add(s1_numerator, 1);
                      int s1 = mul(s1_numerator, mul(inv2, inverse(q2_3)));
                      for (int lam = 0; lam < P; ++lam) {
                        int mu = 0;
                        mu = add(mu, mul(18, mul(mul(mul(q1, q1), q2_2), r0)));
                        mu = add(mu, mul(12, mul(mul(mul(q0, q2_3), r0), 1)));
                        mu = add(mu, mul(12, mul(mul(mul(mul(q1, q1), q1), q2), r1)));
                        mu = add(mu, mul(36, mul(mul(mul(mul(q0, q1), q2_2), r1), 1)));
                        mu = add(mu, mul(3, mul(q2_2, mul(r1, r1))));
                        mu = sub(mu, mul(2, mul(q2_3, s0)));
                        mu = sub(mu, mul(6, mul(mul(mul(q1, q2_2), s1), 1)));
                        mu = sub(mu, mul(2, lam));
                        mu = sub(mu, 2);
                        if (!verify_identity(q0, q1, q2, r0, r1, s0, s1, lam, mu)) continue;
                        ++identities;
                        int generic = lam != 0 && lam != 1 && mu != 0 && mu != 1 && mu != lam;
                        if (generic) ++generic_identities;
                        if ((int)identities <= max_print) {
                            printf(
                                "R3POINT|branch=smooth|q0=%d|q1=%d|q2=%d|u=%d|r0=%d|r1=%d|s0=%d|s1=%d|lambda=%d|mu=%d|generic=%d\n",
                                q0,q1,q2,u,r0,r1,s0,s1,lam,mu,generic
                            );
                        }
                      }
                    }
                }
            }
        }
    }

    uint64_t cusp_outer = 0, cusp_identities = 0, cusp_generic = 0;
    if (run_cusp) for (int q0 = 0; q0 < P; ++q0) {
        for (int q1 = 0; q1 < P; ++q1) {
            for (int q2 = 1; q2 < P; ++q2) {
                int q2_2 = mul(q2, q2), q2_3 = mul(q2_2, q2);
                int r0 = sub(0, mul(q0, q0));
                int s0 = sub(0, mul(mul(q0, q0), q0));
                for (int r1 = 0; r1 < P; ++r1) {
                    ++cusp_outer;
                    int s1_numerator = 0;
                    s1_numerator = add(s1_numerator, mul(12, mul(mul(q1, q2_3), r0)));
                    s1_numerator = add(s1_numerator, mul(18, mul(mul(mul(mul(q1, q1), q2_2), r1), 1)));
                    s1_numerator = add(s1_numerator, mul(12, mul(mul(mul(q0, q2_3), r1), 1)));
                    s1_numerator = add(s1_numerator, 1);
                    int s1 = mul(s1_numerator, mul(inv2, inverse(q2_3)));
                    for (int lam = 0; lam < P; ++lam) {
                        int mu = 0;
                        mu = add(mu, mul(18, mul(mul(mul(q1, q1), q2_2), r0)));
                        mu = add(mu, mul(12, mul(mul(q0, q2_3), r0)));
                        mu = add(mu, mul(12, mul(mul(mul(mul(q1, q1), q1), q2), r1)));
                        mu = add(mu, mul(36, mul(mul(mul(mul(q0, q1), q2_2), r1), 1)));
                        mu = add(mu, mul(3, mul(q2_2, mul(r1, r1))));
                        mu = sub(mu, mul(2, mul(q2_3, s0)));
                        mu = sub(mu, mul(6, mul(mul(mul(q1, q2_2), s1), 1)));
                        mu = sub(mu, mul(2, lam));
                        mu = sub(mu, 2);
                        if (!verify_identity(q0, q1, q2, r0, r1, s0, s1, lam, mu)) continue;
                        ++cusp_identities;
                        int generic = lam != 0 && lam != 1 && mu != 0 && mu != 1 && mu != lam;
                        if (generic) ++cusp_generic;
                        if ((int)(identities + cusp_identities) <= max_print) {
                            printf(
                                "R3POINT|branch=cusp|q0=%d|q1=%d|q2=%d|u=0|r0=%d|r1=%d|s0=%d|s1=%d|lambda=%d|mu=%d|generic=%d\n",
                                q0,q1,q2,r0,r1,s0,s1,lam,mu,generic
                            );
                        }
                    }
                }
            }
        }
    }
    printf(
        "R3SUMMARY|p=%d|smooth_outer=%llu|denominator_zero=%llu|compatible_zero=%llu|smooth_identities=%llu|smooth_generic=%llu|cusp_outer=%llu|cusp_identities=%llu|cusp_generic=%llu|max_print=%d\n",
        P,
        (unsigned long long)outer,
        (unsigned long long)denominator_zero,
        (unsigned long long)compatible_zero,
        (unsigned long long)identities,
        (unsigned long long)generic_identities,
        (unsigned long long)cusp_outer,
        (unsigned long long)cusp_identities,
        (unsigned long long)cusp_generic,
        max_print
    );
    return 0;
}
