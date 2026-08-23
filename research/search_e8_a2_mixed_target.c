/*
 * Exhaustive target-Gram search in a split non-isotrivial E8+A2^3 family.
 *
 * Work over F_p, p > 3.  Put
 *
 *   T = t(t-1),  z = t-lambda,  L = lambda(lambda-1),
 *   k = 3L,       M = kT,
 *
 * and search the elliptic K3 family
 *
 *   y^2 = X^2(X+3M) + c T^2 z^3,              (1)
 *   c = -rho^2/lambda.
 *
 * The residual cubic is R = c z^3 + 4 k^3 T.  We retain exactly the
 * parameters for which R is squarefree.  The resulting fibre configuration
 * is II* + IV@0 + IV@1 + I3@lambda + 3 I1.  The choice k=3L makes the I3
 * tangent lines rational: y = +/- kX.
 *
 * Three section charts encode the target component profiles:
 *
 *   P1: X=zt U1, y=zt V1, deg(U1)<=2, deg(V1)<=4;
 *       nonidentity at 0 and lambda, with V1(lambda)=+kU1(lambda).
 *   P2: X=z U2, y=z V2, deg(U2)<=3, deg(V2)<=5;
 *       nonidentity only at lambda, with V2(lambda)=-kU2(lambda).
 *   P3: deg(X3)<=4, deg(Y3)<=6;
 *       identity at all three finite A2 fibres.
 *
 * The equations after removing the forced factors are
 *
 *   V1^2 = t U1^2 (zU1+3k(t-1)) + c(t-1)^2 z,
 *   V2^2 = U2^2 (zU2+3kT)       + cT^2 z,
 *   Y3^2 = X3^2(X3+3kT)         + cT^2 z^3.
 *
 * We exhaust U1, U2 and X3, reconstruct a polynomial square root, and verify
 * (1) coefficient by coefficient.  A surface is a profile pre-candidate if
 * all three charts are nonempty.
 *
 * For a triple, resolved finite intersections are computed as
 *
 *   gcd(tU1-U2,       tV1-V2),       [remove the forced I3 factor z]
 *   gcd(ztU1-X3,      ztV1-Y3),
 *   gcd(zU2-X3,       zV2-Y3).
 *
 * Requiring all three gcd degrees to be two gives the desired affine
 * intersection numbers.  To certify that no intersection is hidden at the
 * II* fibre at infinity, the program also requires nonzero limiting points
 * (X_4,Y_6) and pairwise distinct limiting points.  In the scaled coordinates
 * (u^4 X,u^6 y), u=1/t, these are distinct smooth points of y^2=X^3 at u=0.
 * Thus the complete triples have P_i.O=0, pairwise intersections 2, and Gram
 *
 *             1/3 [ 8 -1  0 ]
 *                 [ -1 10  0 ].
 *                 [ 0   0 12 ]
 *
 * Output is deterministic JSON Lines.  Counts are exhaustive; examples are
 * print-limited.  A 64-bit FNV-1a audit digest covers every retained section
 * and every affine/complete triple, including records not printed.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wconversion -Wshadow -pedantic \
 *      -o /tmp/search_e8_a2_mixed_target \
 *      research/search_e8_a2_mixed_target.c
 *
 * Run:
 *   /tmp/search_e8_a2_mixed_target 5
 *   /tmp/search_e8_a2_mixed_target 7 20   # at most 20 example records
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define POLY_WIDTH 13
#define P1_U_WIDTH 3
#define P1_V_WIDTH 5
#define P2_U_WIDTH 4
#define P2_V_WIDTH 6
#define X_WIDTH 5
#define Y_WIDTH 7

typedef struct {
    int u[P1_U_WIDTH];
    int v[P1_V_WIDTH];
    int rho;
} P1Section;

typedef struct {
    int u[P2_U_WIDTH];
    int v[P2_V_WIDTH];
} P2Section;

typedef struct {
    int x[X_WIDTH];
    int y[Y_WIDTH];
} P3Section;

typedef struct {
    uint64_t parameter_surfaces;
    uint64_t squarefree_surfaces;
    uint64_t p1_tests;
    uint64_t p2_tests;
    uint64_t p3_tests;
    uint64_t p1_sections;
    uint64_t p2_sections;
    uint64_t p3_sections;
    uint64_t pre_candidate_surfaces;
    uint64_t p1_p2_intersection_two_pairs;
    uint64_t affine_gram_triples;
    uint64_t complete_gram_triples;
    uint64_t complete_surfaces;
} Counts;

static int prime;
static uint64_t audit_hash = UINT64_C(1469598103934665603);

static int mod_i64(int64_t value) {
    int answer = (int)(value % (int64_t)prime);
    return answer < 0 ? answer + prime : answer;
}

static int add_mod(int a, int b) {
    return mod_i64((int64_t)a + (int64_t)b);
}

static int sub_mod(int a, int b) {
    return mod_i64((int64_t)a - (int64_t)b);
}

static int mul_mod(int a, int b) {
    return mod_i64((int64_t)a * (int64_t)b);
}

static int neg_mod(int a) { return a == 0 ? 0 : prime - a; }

static int power_mod(int base, int exponent) {
    int answer = 1;
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            answer = mul_mod(answer, base);
        }
        base = mul_mod(base, base);
        exponent >>= 1;
    }
    return answer;
}

static int inverse_mod(int value) {
    if (value == 0) {
        fputs("internal error: inversion of zero\n", stderr);
        exit(3);
    }
    return power_mod(value, prime - 2);
}

static int is_prime_integer(int value) {
    int divisor;
    if (value < 2) {
        return 0;
    }
    if ((value & 1) == 0) {
        return value == 2;
    }
    for (divisor = 3;
         (int64_t)divisor * (int64_t)divisor <= (int64_t)value;
         divisor += 2) {
        if ((value % divisor) == 0) {
            return 0;
        }
    }
    return 1;
}

static uint64_t integer_power_u64(int base, int exponent) {
    uint64_t result = UINT64_C(1);
    int i;
    for (i = 0; i < exponent; ++i) {
        if (result > UINT64_MAX / (uint64_t)base) {
            fputs("search-space size overflows uint64_t\n", stderr);
            exit(2);
        }
        result *= (uint64_t)base;
    }
    return result;
}

static void audit_integer(int value) {
    uint32_t encoded = (uint32_t)value;
    int shift;
    for (shift = 0; shift < 32; shift += 8) {
        audit_hash ^= (uint64_t)((encoded >> shift) & UINT32_C(255));
        audit_hash *= UINT64_C(1099511628211);
    }
}

static void audit_tag(int tag) { audit_integer(tag); }

static void audit_polynomial(const int *a, int width) {
    int i;
    audit_integer(width);
    for (i = 0; i < width; ++i) {
        audit_integer(a[i]);
    }
}

static void poly_zero(int *a) {
    memset(a, 0, (size_t)POLY_WIDTH * sizeof(*a));
}

static int poly_degree(const int *a) {
    int i;
    for (i = POLY_WIDTH - 1; i >= 0; --i) {
        if (a[i] != 0) {
            return i;
        }
    }
    return -1;
}

static int poly_evaluate(const int *a, int width, int value) {
    int i;
    int answer = 0;
    for (i = width - 1; i >= 0; --i) {
        answer = add_mod(mul_mod(answer, value), a[i]);
    }
    return answer;
}

static void poly_add_scaled(int *target, const int *source, int scalar) {
    int i;
    for (i = 0; i < POLY_WIDTH; ++i) {
        target[i] = add_mod(target[i], mul_mod(scalar, source[i]));
    }
}

static void poly_multiply(const int *a, const int *b, int *answer) {
    int i;
    int j;
    poly_zero(answer);
    for (i = 0; i < POLY_WIDTH; ++i) {
        if (a[i] == 0) {
            continue;
        }
        for (j = 0; i + j < POLY_WIDTH; ++j) {
            if (b[j] != 0) {
                answer[i + j] = add_mod(answer[i + j],
                                        mul_mod(a[i], b[j]));
            }
        }
    }
}

static void poly_derivative(const int *a, int *answer) {
    int i;
    poly_zero(answer);
    for (i = 1; i < POLY_WIDTH; ++i) {
        answer[i - 1] = mul_mod(i, a[i]);
    }
}

static void poly_copy_width(const int *source, int width, int *target) {
    int i;
    poly_zero(target);
    for (i = 0; i < width; ++i) {
        target[i] = source[i];
    }
}

static int polynomial_gcd_degree(const int *left, const int *right) {
    int a[POLY_WIDTH];
    int b[POLY_WIDTH];
    int remainder[POLY_WIDTH];
    int degree_a;
    int degree_b;

    memcpy(a, left, sizeof(a));
    memcpy(b, right, sizeof(b));
    degree_a = poly_degree(a);
    degree_b = poly_degree(b);
    if (degree_a < 0) {
        return degree_b;
    }
    if (degree_b < 0) {
        return degree_a;
    }

    while (degree_b >= 0) {
        int work[POLY_WIDTH];
        int inverse_leading = inverse_mod(b[degree_b]);
        int degree_work;
        memcpy(work, a, sizeof(work));
        degree_work = degree_a;
        while (degree_work >= degree_b) {
            int factor = mul_mod(work[degree_work], inverse_leading);
            int offset = degree_work - degree_b;
            int j;
            for (j = 0; j <= degree_b; ++j) {
                work[offset + j] = sub_mod(
                    work[offset + j], mul_mod(factor, b[j]));
            }
            degree_work = poly_degree(work);
        }
        memcpy(remainder, work, sizeof(remainder));
        memcpy(a, b, sizeof(a));
        memcpy(b, remainder, sizeof(b));
        degree_a = degree_b;
        degree_b = poly_degree(b);
    }
    return degree_a;
}

/* Reconstruct one of the two polynomial square roots, then verify it. */
static int polynomial_square_root(const int *square, int max_degree,
                                  int *root) {
    int degree = poly_degree(square);
    int half_degree;
    int leading_root = -1;
    int candidate;
    int j;
    int coefficient;

    poly_zero(root);
    if (degree < 0) {
        return 1;
    }
    if ((degree & 1) != 0) {
        return 0;
    }
    half_degree = degree / 2;
    if (half_degree > max_degree) {
        return 0;
    }
    for (candidate = 1; candidate < prime; ++candidate) {
        if (mul_mod(candidate, candidate) == square[degree]) {
            leading_root = candidate;
            break;
        }
    }
    if (leading_root < 0) {
        return 0;
    }
    root[half_degree] = leading_root;

    for (j = half_degree - 1; j >= 0; --j) {
        int known = 0;
        int target_degree = half_degree + j;
        int i;
        for (i = j + 1; i <= half_degree; ++i) {
            int other = target_degree - i;
            if (other >= 0 && other <= half_degree && other != j) {
                known = add_mod(known, mul_mod(root[i], root[other]));
            }
        }
        root[j] = mul_mod(sub_mod(square[target_degree], known),
                          inverse_mod(mul_mod(2, root[half_degree])));
    }

    for (coefficient = 0; coefficient < POLY_WIDTH; ++coefficient) {
        int observed = 0;
        int i;
        for (i = 0; i <= max_degree; ++i) {
            int other = coefficient - i;
            if (other >= 0 && other <= max_degree) {
                observed = add_mod(observed, mul_mod(root[i], root[other]));
            }
        }
        if (observed != square[coefficient]) {
            return 0;
        }
    }
    return 1;
}

static void decode_polynomial(uint64_t code, int width, int *answer) {
    int i;
    for (i = 0; i < width; ++i) {
        answer[i] = (int)(code % (uint64_t)prime);
        code /= (uint64_t)prime;
    }
}

static void negate_polynomial(int *a, int width) {
    int i;
    for (i = 0; i < width; ++i) {
        a[i] = neg_mod(a[i]);
    }
}

static void build_basic_polynomials(int lambda, int *t_poly, int *tm1,
                                    int *z, int *T, int *T2, int *z2,
                                    int *z3, int *zt) {
    int scratch[POLY_WIDTH];
    poly_zero(t_poly);
    poly_zero(tm1);
    poly_zero(z);
    t_poly[1] = 1;
    tm1[0] = neg_mod(1);
    tm1[1] = 1;
    z[0] = neg_mod(lambda);
    z[1] = 1;
    poly_multiply(t_poly, tm1, T);
    poly_multiply(T, T, T2);
    poly_multiply(z, z, z2);
    poly_multiply(z2, z, z3);
    poly_multiply(z, t_poly, scratch);
    memcpy(zt, scratch, sizeof(scratch));
}

static void build_p1_square(const int *u_small, int c, int k,
                            const int *t_poly, const int *tm1, const int *z,
                            int *answer) {
    int u[POLY_WIDTH];
    int u2[POLY_WIDTH];
    int zu[POLY_WIDTH];
    int inner[POLY_WIDTH];
    int product[POLY_WIDTH];
    int scratch[POLY_WIDTH];
    int tm1_squared[POLY_WIDTH];
    int correction[POLY_WIDTH];

    poly_copy_width(u_small, P1_U_WIDTH, u);
    poly_multiply(u, u, u2);
    poly_multiply(z, u, zu);
    memcpy(inner, zu, sizeof(inner));
    poly_add_scaled(inner, tm1, mul_mod(3, k));
    poly_multiply(u2, inner, product);
    poly_multiply(t_poly, product, scratch);
    poly_multiply(tm1, tm1, tm1_squared);
    poly_multiply(tm1_squared, z, correction);
    poly_add_scaled(scratch, correction, c);
    memcpy(answer, scratch, sizeof(scratch));
}

static void build_p2_square(const int *u_small, int c, int k,
                            const int *z, const int *T, const int *T2,
                            int *answer) {
    int u[POLY_WIDTH];
    int u2[POLY_WIDTH];
    int zu[POLY_WIDTH];
    int inner[POLY_WIDTH];
    int product[POLY_WIDTH];
    int correction[POLY_WIDTH];

    poly_copy_width(u_small, P2_U_WIDTH, u);
    poly_multiply(u, u, u2);
    poly_multiply(z, u, zu);
    memcpy(inner, zu, sizeof(inner));
    poly_add_scaled(inner, T, mul_mod(3, k));
    poly_multiply(u2, inner, product);
    poly_multiply(T2, z, correction);
    poly_add_scaled(product, correction, c);
    memcpy(answer, product, sizeof(product));
}

static void build_p3_square(const int *x_small, int c, int k,
                            const int *T, const int *T2, const int *z3,
                            int *answer) {
    int x[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int inner[POLY_WIDTH];
    int product[POLY_WIDTH];
    int correction[POLY_WIDTH];

    poly_copy_width(x_small, X_WIDTH, x);
    poly_multiply(x, x, x2);
    memcpy(inner, x, sizeof(inner));
    poly_add_scaled(inner, T, mul_mod(3, k));
    poly_multiply(x2, inner, product);
    poly_multiply(T2, z3, correction);
    poly_add_scaled(product, correction, c);
    memcpy(answer, product, sizeof(product));
}

static void multiply_small_by_poly(const int *small, int small_width,
                                   const int *factor, int *answer) {
    int expanded[POLY_WIDTH];
    poly_copy_width(small, small_width, expanded);
    poly_multiply(expanded, factor, answer);
}

static int verify_full_section(const int *x_small, int x_width,
                               const int *y_small, int y_width,
                               int c, int k, const int *T, const int *T2,
                               const int *z3) {
    int x[POLY_WIDTH];
    int y[POLY_WIDTH];
    int y2[POLY_WIDTH];
    int expected[POLY_WIDTH];
    poly_copy_width(x_small, x_width, x);
    poly_copy_width(y_small, y_width, y);
    poly_multiply(y, y, y2);
    build_p3_square(x_small, c, k, T, T2, z3, expected);
    return memcmp(y2, expected, sizeof(y2)) == 0;
}

static int residual_is_squarefree(int c, int k, const int *T,
                                  const int *z3) {
    int residual[POLY_WIDTH];
    int derivative[POLY_WIDTH];
    poly_zero(residual);
    poly_add_scaled(residual, z3, c);
    poly_add_scaled(residual, T, mul_mod(4, mul_mod(mul_mod(k, k), k)));
    if (poly_degree(residual) != 3) {
        return 0;
    }
    poly_derivative(residual, derivative);
    return polynomial_gcd_degree(residual, derivative) == 0;
}

static size_t checked_capacity(uint64_t count, size_t element_size) {
    if (count > (uint64_t)(SIZE_MAX / element_size)) {
        fputs("section storage size overflows size_t\n", stderr);
        exit(2);
    }
    return (size_t)count;
}

static void *checked_calloc(size_t count, size_t size) {
    void *answer = calloc(count, size);
    if (answer == NULL && count != 0U) {
        fputs("failed to allocate exhaustive section storage\n", stderr);
        exit(2);
    }
    return answer;
}

static size_t enumerate_p1(int lambda, int c, int k,
                           const int *t_poly, const int *tm1, const int *z,
                           const int *T, const int *T2, const int *z3,
                           const int *zt, P1Section *sections,
                           Counts *counts) {
    uint64_t limit = integer_power_u64(prime, P1_U_WIDTH);
    size_t found = 0U;
    uint64_t code;
    for (code = 0; code < limit; ++code) {
        int u[P1_U_WIDTH];
        int square[POLY_WIDTH];
        int root[POLY_WIDTH];
        int full_x[POLY_WIDTH];
        int full_y[POLY_WIDTH];
        int value_u_lambda;
        int value_v_lambda;
        decode_polynomial(code, P1_U_WIDTH, u);
        ++counts->p1_tests;
        if (poly_evaluate(u, P1_U_WIDTH, 1) == 0) {
            continue;
        }
        value_u_lambda = poly_evaluate(u, P1_U_WIDTH, lambda);
        if (value_u_lambda == 0) {
            continue;
        }
        build_p1_square(u, c, k, t_poly, tm1, z, square);
        if (!polynomial_square_root(square, 4, root)) {
            continue;
        }
        value_v_lambda = poly_evaluate(root, P1_V_WIDTH, lambda);
        if (value_v_lambda != mul_mod(k, value_u_lambda)) {
            negate_polynomial(root, P1_V_WIDTH);
            value_v_lambda = poly_evaluate(root, P1_V_WIDTH, lambda);
        }
        if (value_v_lambda != mul_mod(k, value_u_lambda)) {
            continue;
        }
        if (root[0] == 0 || c != neg_mod(mul_mod(mul_mod(root[0], root[0]),
                                                 inverse_mod(lambda)))) {
            fputs("internal error: P1 rho parametrization failed\n", stderr);
            exit(3);
        }
        multiply_small_by_poly(u, P1_U_WIDTH, zt, full_x);
        multiply_small_by_poly(root, P1_V_WIDTH, zt, full_y);
        if (!verify_full_section(full_x, X_WIDTH, full_y, Y_WIDTH,
                                 c, k, T, T2, z3)) {
            fputs("internal error: P1 failed full section identity\n", stderr);
            exit(3);
        }
        memcpy(sections[found].u, u, sizeof(u));
        memcpy(sections[found].v, root,
               (size_t)P1_V_WIDTH * sizeof(*root));
        sections[found].rho = root[0];
        audit_tag(101);
        audit_integer(lambda);
        audit_integer(c);
        audit_polynomial(u, P1_U_WIDTH);
        audit_polynomial(root, P1_V_WIDTH);
        ++found;
        ++counts->p1_sections;
    }
    return found;
}

static size_t enumerate_p2(int lambda, int c, int k,
                           const int *z, const int *T, const int *T2,
                           const int *z3, P2Section *sections,
                           Counts *counts) {
    uint64_t limit = integer_power_u64(prime, P2_U_WIDTH);
    size_t found = 0U;
    uint64_t code;
    for (code = 0; code < limit; ++code) {
        int u[P2_U_WIDTH];
        int square[POLY_WIDTH];
        int root[POLY_WIDTH];
        int full_x[POLY_WIDTH];
        int full_y[POLY_WIDTH];
        int value_u_lambda;
        int value_v_lambda;
        decode_polynomial(code, P2_U_WIDTH, u);
        ++counts->p2_tests;
        if (u[0] == 0 || poly_evaluate(u, P2_U_WIDTH, 1) == 0) {
            continue;
        }
        value_u_lambda = poly_evaluate(u, P2_U_WIDTH, lambda);
        if (value_u_lambda == 0) {
            continue;
        }
        build_p2_square(u, c, k, z, T, T2, square);
        if (!polynomial_square_root(square, 5, root)) {
            continue;
        }
        value_v_lambda = poly_evaluate(root, P2_V_WIDTH, lambda);
        if (value_v_lambda != neg_mod(mul_mod(k, value_u_lambda))) {
            negate_polynomial(root, P2_V_WIDTH);
            value_v_lambda = poly_evaluate(root, P2_V_WIDTH, lambda);
        }
        if (value_v_lambda != neg_mod(mul_mod(k, value_u_lambda))) {
            continue;
        }
        multiply_small_by_poly(u, P2_U_WIDTH, z, full_x);
        multiply_small_by_poly(root, P2_V_WIDTH, z, full_y);
        if (!verify_full_section(full_x, X_WIDTH, full_y, Y_WIDTH,
                                 c, k, T, T2, z3)) {
            fputs("internal error: P2 failed full section identity\n", stderr);
            exit(3);
        }
        memcpy(sections[found].u, u, sizeof(u));
        memcpy(sections[found].v, root,
               (size_t)P2_V_WIDTH * sizeof(*root));
        audit_tag(102);
        audit_integer(lambda);
        audit_integer(c);
        audit_polynomial(u, P2_U_WIDTH);
        audit_polynomial(root, P2_V_WIDTH);
        ++found;
        ++counts->p2_sections;
    }
    return found;
}

static size_t enumerate_p3(int lambda, int c, int k,
                           const int *T, const int *T2, const int *z3,
                           P3Section *sections, Counts *counts) {
    uint64_t limit = integer_power_u64(prime, X_WIDTH);
    size_t found = 0U;
    uint64_t code;
    for (code = 0; code < limit; ++code) {
        int x[X_WIDTH];
        int square[POLY_WIDTH];
        int root[POLY_WIDTH];
        int sign;
        decode_polynomial(code, X_WIDTH, x);
        ++counts->p3_tests;
        if (x[0] == 0 || poly_evaluate(x, X_WIDTH, 1) == 0
            || poly_evaluate(x, X_WIDTH, lambda) == 0) {
            continue;
        }
        build_p3_square(x, c, k, T, T2, z3, square);
        if (!polynomial_square_root(square, 6, root)) {
            continue;
        }
        for (sign = 0; sign < 2; ++sign) {
            if (sign != 0) {
                negate_polynomial(root, Y_WIDTH);
            }
            if (!verify_full_section(x, X_WIDTH, root, Y_WIDTH,
                                     c, k, T, T2, z3)) {
                fputs("internal error: P3 failed full section identity\n",
                      stderr);
                exit(3);
            }
            memcpy(sections[found].x, x, sizeof(x));
            memcpy(sections[found].y, root,
                   (size_t)Y_WIDTH * sizeof(*root));
            audit_tag(103);
            audit_integer(lambda);
            audit_integer(c);
            audit_polynomial(x, X_WIDTH);
            audit_polynomial(root, Y_WIDTH);
            ++found;
            ++counts->p3_sections;
        }
        /* The second iteration left root negated; no later use is made of it. */
    }
    return found;
}

static int p1_p2_intersection_degree(const P1Section *p1,
                                     const P2Section *p2) {
    int dx[POLY_WIDTH];
    int dy[POLY_WIDTH];
    int i;
    poly_zero(dx);
    poly_zero(dy);
    for (i = 0; i < P1_U_WIDTH; ++i) {
        dx[i + 1] = add_mod(dx[i + 1], p1->u[i]);
    }
    for (i = 0; i < P2_U_WIDTH; ++i) {
        dx[i] = sub_mod(dx[i], p2->u[i]);
    }
    for (i = 0; i < P1_V_WIDTH; ++i) {
        dy[i + 1] = add_mod(dy[i + 1], p1->v[i]);
    }
    for (i = 0; i < P2_V_WIDTH; ++i) {
        dy[i] = sub_mod(dy[i], p2->v[i]);
    }
    return polynomial_gcd_degree(dx, dy);
}

static int p1_p3_intersection_degree(const P1Section *p1,
                                     const P3Section *p3,
                                     const int *zt) {
    int x1[POLY_WIDTH];
    int y1[POLY_WIDTH];
    int dx[POLY_WIDTH];
    int dy[POLY_WIDTH];
    int i;
    multiply_small_by_poly(p1->u, P1_U_WIDTH, zt, x1);
    multiply_small_by_poly(p1->v, P1_V_WIDTH, zt, y1);
    memcpy(dx, x1, sizeof(dx));
    memcpy(dy, y1, sizeof(dy));
    for (i = 0; i < X_WIDTH; ++i) {
        dx[i] = sub_mod(dx[i], p3->x[i]);
    }
    for (i = 0; i < Y_WIDTH; ++i) {
        dy[i] = sub_mod(dy[i], p3->y[i]);
    }
    return polynomial_gcd_degree(dx, dy);
}

static int p2_p3_intersection_degree(const P2Section *p2,
                                     const P3Section *p3,
                                     const int *z) {
    int x2[POLY_WIDTH];
    int y2[POLY_WIDTH];
    int dx[POLY_WIDTH];
    int dy[POLY_WIDTH];
    int i;
    multiply_small_by_poly(p2->u, P2_U_WIDTH, z, x2);
    multiply_small_by_poly(p2->v, P2_V_WIDTH, z, y2);
    memcpy(dx, x2, sizeof(dx));
    memcpy(dy, y2, sizeof(dy));
    for (i = 0; i < X_WIDTH; ++i) {
        dx[i] = sub_mod(dx[i], p3->x[i]);
    }
    for (i = 0; i < Y_WIDTH; ++i) {
        dy[i] = sub_mod(dy[i], p3->y[i]);
    }
    return polynomial_gcd_degree(dx, dy);
}

static void p1_full(const P1Section *section, const int *zt,
                    int *x, int *y) {
    multiply_small_by_poly(section->u, P1_U_WIDTH, zt, x);
    multiply_small_by_poly(section->v, P1_V_WIDTH, zt, y);
}

static void p2_full(const P2Section *section, const int *z,
                    int *x, int *y) {
    multiply_small_by_poly(section->u, P2_U_WIDTH, z, x);
    multiply_small_by_poly(section->v, P2_V_WIDTH, z, y);
}

static int limiting_point_is_smooth(int x4, int y6) {
    return x4 != 0 && y6 != 0 && mul_mod(y6, y6) == mul_mod(x4, mul_mod(x4, x4));
}

static int limiting_points_differ(int x1, int y1, int x2, int y2) {
    return x1 != x2 || y1 != y2;
}

static int infinity_is_safe(const P1Section *p1, const P2Section *p2,
                            const P3Section *p3, const int *zt,
                            const int *z) {
    int x1[POLY_WIDTH];
    int y1[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int y2[POLY_WIDTH];
    int x3 = p3->x[4];
    int y3 = p3->y[6];
    p1_full(p1, zt, x1, y1);
    p2_full(p2, z, x2, y2);
    if (!limiting_point_is_smooth(x1[4], y1[6])
        || !limiting_point_is_smooth(x2[4], y2[6])
        || !limiting_point_is_smooth(x3, y3)) {
        return 0;
    }
    return limiting_points_differ(x1[4], y1[6], x2[4], y2[6])
           && limiting_points_differ(x1[4], y1[6], x3, y3)
           && limiting_points_differ(x2[4], y2[6], x3, y3);
}

static void print_array(const int *values, int width) {
    int i;
    putchar('[');
    for (i = 0; i < width; ++i) {
        printf("%s%d", i == 0 ? "" : ",", values[i]);
    }
    putchar(']');
}

static void print_chart_surface(int lambda, int c, int k,
                                const P1Section *p1, size_t p1_count,
                                const P2Section *p2, size_t p2_count,
                                const P3Section *p3, size_t p3_count) {
    printf("{\"event\":\"nonempty_chart_surface\",\"prime\":%d,"
           "\"lambda\":%d,\"c\":%d,\"k\":%d,"
           "\"section_counts\":[%zu,%zu,%zu],\"first\":{",
           prime, lambda, c, k, p1_count, p2_count, p3_count);
    if (p1_count != 0U) {
        printf("\"P1\":{\"rho\":%d,\"U\":", p1[0].rho);
        print_array(p1[0].u, P1_U_WIDTH);
        printf(",\"V\":");
        print_array(p1[0].v, P1_V_WIDTH);
        putchar('}');
    } else {
        printf("\"P1\":null");
    }
    if (p2_count != 0U) {
        printf(",\"P2\":{\"U\":");
        print_array(p2[0].u, P2_U_WIDTH);
        printf(",\"V\":");
        print_array(p2[0].v, P2_V_WIDTH);
        putchar('}');
    } else {
        printf(",\"P2\":null");
    }
    if (p3_count != 0U) {
        printf(",\"P3\":{\"X\":");
        print_array(p3[0].x, X_WIDTH);
        printf(",\"Y\":");
        print_array(p3[0].y, Y_WIDTH);
        putchar('}');
    } else {
        printf(",\"P3\":null");
    }
    printf("}}\n");
}

static void print_complete_example(int lambda, int c, int k,
                                   const P1Section *p1,
                                   const P2Section *p2,
                                   const P3Section *p3) {
    printf("{\"event\":\"complete_gram_example\",\"prime\":%d,"
           "\"lambda\":%d,\"c\":%d,\"k\":%d,\"rho\":%d,"
           "\"P1\":{\"U\":", prime, lambda, c, k, p1->rho);
    print_array(p1->u, P1_U_WIDTH);
    printf(",\"V\":");
    print_array(p1->v, P1_V_WIDTH);
    printf("},\"P2\":{\"U\":");
    print_array(p2->u, P2_U_WIDTH);
    printf(",\"V\":");
    print_array(p2->v, P2_V_WIDTH);
    printf("},\"P3\":{\"X\":");
    print_array(p3->x, X_WIDTH);
    printf(",\"Y\":");
    print_array(p3->y, Y_WIDTH);
    printf("},\"resolved_intersections\":[2,2,2],"
           "\"infinity_intersections\":[0,0,0],"
           "\"gram_numerator_over_3\":[[8,-1,0],[-1,10,0],[0,0,12]]}\n");
}

static void audit_triple(int tag, int lambda, int c, size_t i, size_t j,
                         size_t m) {
    audit_tag(tag);
    audit_integer(lambda);
    audit_integer(c);
    audit_integer((int)i);
    audit_integer((int)j);
    audit_integer((int)m);
}

int main(int argc, char **argv) {
    Counts counts = {0};
    int max_print = 20;
    int printed = 0;
    uint64_t p3_capacity_u64;
    size_t p1_capacity;
    size_t p2_capacity;
    size_t p3_capacity;
    P1Section *p1_sections;
    P2Section *p2_sections;
    P3Section *p3_sections;
    int lambda;

    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: %s PRIME [MAX_PRINT]\n", argv[0]);
        return 2;
    }
    prime = atoi(argv[1]);
    if (!is_prime_integer(prime) || prime <= 3) {
        fputs("PRIME must be a prime integer greater than 3\n", stderr);
        return 2;
    }
    if (prime > 31) {
        fputs("PRIME above 31 is refused: the exhaustive p^5 storage/search "
              "is intentionally bounded\n", stderr);
        return 2;
    }
    if (argc == 3) {
        max_print = atoi(argv[2]);
        if (max_print < 0) {
            fputs("MAX_PRINT must be nonnegative\n", stderr);
            return 2;
        }
    }

    p1_capacity = checked_capacity(integer_power_u64(prime, P1_U_WIDTH),
                                   sizeof(*p1_sections));
    p2_capacity = checked_capacity(integer_power_u64(prime, P2_U_WIDTH),
                                   sizeof(*p2_sections));
    p3_capacity_u64 = 2U * integer_power_u64(prime, X_WIDTH);
    p3_capacity = checked_capacity(p3_capacity_u64, sizeof(*p3_sections));
    p1_sections = checked_calloc(p1_capacity, sizeof(*p1_sections));
    p2_sections = checked_calloc(p2_capacity, sizeof(*p2_sections));
    p3_sections = checked_calloc(p3_capacity, sizeof(*p3_sections));

    printf("{\"event\":\"start\",\"prime\":%d,\"max_print\":%d,"
           "\"family\":\"y^2=X^2(X+3kT)+cT^2z^3\","
           "\"k_rule\":\"k=3lambda(lambda-1)\","
           "\"c_rule\":\"-c*lambda is a nonzero square\"}\n",
           prime, max_print);

    for (lambda = 2; lambda < prime; ++lambda) {
        int t_poly[POLY_WIDTH];
        int tm1[POLY_WIDTH];
        int z[POLY_WIDTH];
        int T[POLY_WIDTH];
        int T2[POLY_WIDTH];
        int z2[POLY_WIDTH];
        int z3[POLY_WIDTH];
        int zt[POLY_WIDTH];
        int L = mul_mod(lambda, sub_mod(lambda, 1));
        int k = mul_mod(3, L);
        int c;
        build_basic_polynomials(lambda, t_poly, tm1, z, T, T2, z2, z3,
                                zt);
        (void)z2;

        for (c = 1; c < prime; ++c) {
            int rho_square = neg_mod(mul_mod(c, lambda));
            size_t p1_count;
            size_t p2_count;
            size_t p3_count;
            uint64_t complete_before;
            size_t i;
            if (power_mod(rho_square, (prime - 1) / 2) != 1) {
                continue;
            }
            ++counts.parameter_surfaces;
            if (!residual_is_squarefree(c, k, T, z3)) {
                continue;
            }
            ++counts.squarefree_surfaces;
            audit_tag(100);
            audit_integer(lambda);
            audit_integer(c);
            audit_integer(k);

            p1_count = enumerate_p1(lambda, c, k, t_poly, tm1, z, T, T2,
                                    z3, zt, p1_sections, &counts);
            p2_count = enumerate_p2(lambda, c, k, z, T, T2, z3,
                                    p2_sections, &counts);
            p3_count = enumerate_p3(lambda, c, k, T, T2, z3,
                                    p3_sections, &counts);
            if (p1_count != 0U || p2_count != 0U || p3_count != 0U) {
                print_chart_surface(lambda, c, k, p1_sections, p1_count,
                                    p2_sections, p2_count,
                                    p3_sections, p3_count);
            }
            if (p1_count == 0U || p2_count == 0U || p3_count == 0U) {
                continue;
            }
            ++counts.pre_candidate_surfaces;
            complete_before = counts.complete_gram_triples;

            for (i = 0U; i < p1_count; ++i) {
                size_t j;
                for (j = 0U; j < p2_count; ++j) {
                    size_t m;
                    if (p1_p2_intersection_degree(&p1_sections[i],
                                                  &p2_sections[j]) != 2) {
                        continue;
                    }
                    ++counts.p1_p2_intersection_two_pairs;
                    for (m = 0U; m < p3_count; ++m) {
                        if (p1_p3_intersection_degree(&p1_sections[i],
                                                      &p3_sections[m], zt)
                                != 2
                            || p2_p3_intersection_degree(&p2_sections[j],
                                                         &p3_sections[m], z)
                                   != 2) {
                            continue;
                        }
                        ++counts.affine_gram_triples;
                        audit_triple(201, lambda, c, i, j, m);
                        if (!infinity_is_safe(&p1_sections[i],
                                              &p2_sections[j],
                                              &p3_sections[m], zt, z)) {
                            continue;
                        }
                        ++counts.complete_gram_triples;
                        audit_triple(202, lambda, c, i, j, m);
                        if (printed < max_print) {
                            print_complete_example(lambda, c, k,
                                                   &p1_sections[i],
                                                   &p2_sections[j],
                                                   &p3_sections[m]);
                            ++printed;
                        }
                    }
                }
            }
            if (counts.complete_gram_triples > complete_before) {
                ++counts.complete_surfaces;
            }
            printf("{\"event\":\"pre_candidate_surface\",\"prime\":%d,"
                   "\"lambda\":%d,\"c\":%d,\"k\":%d,"
                   "\"section_counts\":[%zu,%zu,%zu],"
                   "\"complete_triples_on_surface\":%" PRIu64 "}\n",
                   prime, lambda, c, k, p1_count, p2_count, p3_count,
                   counts.complete_gram_triples - complete_before);
        }
    }

    printf("{\"event\":\"summary\",\"prime\":%d,"
           "\"parameter_surfaces\":%" PRIu64 ","
           "\"squarefree_surfaces\":%" PRIu64 ","
           "\"tests\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"sections\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"pre_candidate_surfaces\":%" PRIu64 ","
           "\"P1_P2_intersection_two_pairs\":%" PRIu64 ","
           "\"affine_gram_triples\":%" PRIu64 ","
           "\"complete_gram_triples\":%" PRIu64 ","
           "\"complete_surfaces\":%" PRIu64 ","
           "\"printed_complete_examples\":%d,"
           "\"audit_fnv1a64\":\"%016" PRIx64 "\","
           "\"claim_boundary\":\"complete triples certify the target Gram "
           "over F_p; lifting to characteristic zero remains separate\"}\n",
           prime, counts.parameter_surfaces, counts.squarefree_surfaces,
           counts.p1_tests, counts.p2_tests, counts.p3_tests,
           counts.p1_sections, counts.p2_sections, counts.p3_sections,
           counts.pre_candidate_surfaces,
           counts.p1_p2_intersection_two_pairs,
           counts.affine_gram_triples, counts.complete_gram_triples,
           counts.complete_surfaces, printed, audit_hash);

    free(p1_sections);
    free(p2_sections);
    free(p3_sections);
    return 0;
}
