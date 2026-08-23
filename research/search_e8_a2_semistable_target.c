/*
 * Exhaustive target-Gram search on a rational split cover of the full
 * II* + 3 I3 semistable E8+A2^3 chart.
 *
 * Work over F_p, p>3.  Let D=t(t-1)(t-lambda).  The shifted model is
 *
 *   y^2 = X^2(X+3a) - 6D beta X + D^2 gamma,                 (1)
 *
 * where deg(a,beta,gamma)<=(2,1,1) and
 *
 *   a gamma - 3 beta^2 = dD.
 *
 * The split cover is parametrized by (lambda,m,R,W), after using the
 * Weierstrass homothety to set sqrt(3a(0))=1.  Put
 *
 *   x = (m^2-2m+lambda)/(m^2-lambda),
 *   y = (-m^2+2m lambda-lambda)/(m^2-lambda),
 *   q = ((1-lambda)+lambda R x)/y.
 *
 * Then y^2=(1-lambda)+lambda*x^2, and interpolation from
 *
 *   3a(0)=1,       3a(1)=R^2,       3a(lambda)=q^2,
 *   3beta(0)=W,    3beta(1)=R*x*W,
 *   gamma(0)=W^2, gamma(1)=x^2*W^2
 *
 * gives the required polynomial identity.  The program checks that identity
 * coefficient-by-coefficient, rather than trusting the parametrization.
 *
 * We retain only the exact Kodaira open set: gamma_1!=0 and
 *
 *   H=4a^2d+12a beta gamma-32beta^3+D gamma^2
 *
 * is degree five, squarefree and coprime to D.  Thus
 * Delta=-432 D^3 H has II*+3I3+5I1.
 *
 * Section charts and target profiles are
 *
 *   P1: X=t(t-lambda)U1, y=t(t-lambda)V1, deg(U1,V1)<=(2,4),
 *       nonidentity at 0 and lambda, positive q-branch at lambda;
 *   P2: X=(t-lambda)U2, y=(t-lambda)V2, deg(U2,V2)<=(3,5),
 *       nonidentity only at lambda, negative q-branch there;
 *   P3: deg(X3,Y3)<=(4,6), identity at all three I3 fibres.
 *
 * The opposite lambda branches give local pairing 1/3.  Ambiguous deeper
 * nodal contacts are deliberately saturated away.  The cascade exhausts P1
 * first, P2 only on P1-nonempty surfaces, and P3 only after a P1/P2 pair has
 * resolved intersection degree two.
 *
 * Finite resolved intersections use
 *
 *   gcd(tU1-U2,tV1-V2),
 *   gcd(t(t-lambda)U1-X3,t(t-lambda)V1-Y3),
 *   gcd((t-lambda)U2-X3,(t-lambda)V2-Y3).
 *
 * Complete triples additionally have three distinct nonzero smooth limiting
 * points on y^2=X^3 at infinity.  They therefore certify P_i.O=0, all three
 * pairwise intersections equal two, and Gram
 *
 *          1/3 [[8,-1,0],[-1,10,0],[0,0,12]]
 *
 * over F_p.  This is finite-field evidence only; characteristic-zero lifting
 * and rank 31 remain separate.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wconversion -Wshadow -pedantic -Werror \
 *      -o /tmp/search_e8_a2_semistable_target \
 *      research/search_e8_a2_semistable_target.c
 *
 * Run:
 *   /tmp/search_e8_a2_semistable_target 5
 *   /tmp/search_e8_a2_semistable_target 7 20
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
    uint64_t raw_parameter_tuples;
    uint64_t rational_chart_tuples;
    uint64_t kodaira_surfaces;
    uint64_t p1_stage_surfaces;
    uint64_t p2_stage_surfaces;
    uint64_t p3_stage_surfaces;
    uint64_t p1_tests;
    uint64_t p2_tests;
    uint64_t p3_tests;
    uint64_t p1_sections;
    uint64_t p2_sections;
    uint64_t p3_sections;
    uint64_t p1_p2_intersection_two_pairs;
    uint64_t pre_candidate_surfaces;
    uint64_t affine_gram_triples;
    uint64_t complete_gram_triples;
    uint64_t complete_surfaces;
} Counts;

typedef struct {
    int lambda;
    int m;
    int R;
    int W;
    int conic_x;
    int conic_y;
    int q;
    int s0;
    int s_lambda;
    int d;
    int a[POLY_WIDTH];
    int beta[POLY_WIDTH];
    int gamma[POLY_WIDTH];
    int D[POLY_WIDTH];
    int D2[POLY_WIDTH];
    int t[POLY_WIDTH];
    int tm1[POLY_WIDTH];
    int z[POLY_WIDTH];
    int tz[POLY_WIDTH];
    int T[POLY_WIDTH];
} Surface;

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

static void audit_u64(uint64_t value) {
    int shift;
    for (shift = 0; shift < 64; shift += 8) {
        audit_hash ^= (value >> shift) & UINT64_C(255);
        audit_hash *= UINT64_C(1099511628211);
    }
}

static void audit_integer(int value) {
    audit_u64((uint64_t)(uint32_t)value);
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
        memcpy(a, b, sizeof(a));
        memcpy(b, work, sizeof(b));
        degree_a = degree_b;
        degree_b = poly_degree(b);
    }
    return degree_a;
}

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

static void multiply_small_by_poly(const int *small, int small_width,
                                   const int *factor, int *answer) {
    int expanded[POLY_WIDTH];
    poly_copy_width(small, small_width, expanded);
    poly_multiply(expanded, factor, answer);
}

static void build_basic_polynomials(Surface *surface) {
    poly_zero(surface->t);
    poly_zero(surface->tm1);
    poly_zero(surface->z);
    surface->t[1] = 1;
    surface->tm1[0] = neg_mod(1);
    surface->tm1[1] = 1;
    surface->z[0] = neg_mod(surface->lambda);
    surface->z[1] = 1;
    poly_multiply(surface->t, surface->tm1, surface->T);
    poly_multiply(surface->T, surface->z, surface->D);
    poly_multiply(surface->D, surface->D, surface->D2);
    poly_multiply(surface->t, surface->z, surface->tz);
}

static int relation_holds(const Surface *surface) {
    int left[POLY_WIDTH];
    int beta2[POLY_WIDTH];
    int right[POLY_WIDTH];
    poly_multiply(surface->a, surface->gamma, left);
    poly_multiply(surface->beta, surface->beta, beta2);
    poly_add_scaled(left, beta2, neg_mod(3));
    poly_zero(right);
    poly_add_scaled(right, surface->D, surface->d);
    return memcmp(left, right, sizeof(left)) == 0;
}

static void build_H(const Surface *surface, int *H) {
    int a2[POLY_WIDTH];
    int ab[POLY_WIDTH];
    int abg[POLY_WIDTH];
    int beta2[POLY_WIDTH];
    int beta3[POLY_WIDTH];
    int gamma2[POLY_WIDTH];
    int Dgamma2[POLY_WIDTH];
    poly_multiply(surface->a, surface->a, a2);
    poly_multiply(surface->a, surface->beta, ab);
    poly_multiply(ab, surface->gamma, abg);
    poly_multiply(surface->beta, surface->beta, beta2);
    poly_multiply(beta2, surface->beta, beta3);
    poly_multiply(surface->gamma, surface->gamma, gamma2);
    poly_multiply(surface->D, gamma2, Dgamma2);
    poly_zero(H);
    poly_add_scaled(H, a2, mul_mod(4, surface->d));
    poly_add_scaled(H, abg, 12);
    poly_add_scaled(H, beta3, neg_mod(32));
    poly_add_scaled(H, Dgamma2, 1);
}

/* Return 0 off the rational chart, 1 on the chart but off the Kodaira open,
 * and 2 on the exact II*+3I3+5I1 open set. */
static int build_surface(int lambda, int m, int R, int W,
                         Surface *surface) {
    int denominator = sub_mod(mul_mod(m, m), lambda);
    int numerator_x;
    int numerator_y;
    int inv_denominator;
    int inv3 = inverse_mod(3);
    int one_minus_lambda = sub_mod(1, lambda);
    int a_at_0 = inv3;
    int a_at_1;
    int a_at_lambda;
    int delta_1;
    int delta_lambda;
    int interpolation_denominator;
    int H[POLY_WIDTH];
    int derivative[POLY_WIDTH];
    if (denominator == 0) {
        return 0;
    }
    inv_denominator = inverse_mod(denominator);
    numerator_x = add_mod(sub_mod(mul_mod(m, m), mul_mod(2, m)), lambda);
    numerator_y = sub_mod(add_mod(neg_mod(mul_mod(m, m)),
                                  mul_mod(mul_mod(2, m), lambda)), lambda);
    memset(surface, 0, sizeof(*surface));
    surface->lambda = lambda;
    surface->m = m;
    surface->R = R;
    surface->W = W;
    surface->conic_x = mul_mod(numerator_x, inv_denominator);
    surface->conic_y = mul_mod(numerator_y, inv_denominator);
    if (surface->conic_y == 0
        || mul_mod(surface->conic_x, surface->conic_x) == 1) {
        return 0;
    }
    if (mul_mod(surface->conic_y, surface->conic_y)
        != add_mod(one_minus_lambda,
                   mul_mod(lambda,
                           mul_mod(surface->conic_x, surface->conic_x)))) {
        fputs("internal error: conic parametrization failed\n", stderr);
        exit(3);
    }
    surface->q = mul_mod(
        add_mod(one_minus_lambda,
                mul_mod(lambda, mul_mod(R, surface->conic_x))),
        inverse_mod(surface->conic_y));
    if (surface->q == 0) {
        return 0;
    }
    surface->s0 = W;
    surface->s_lambda = mul_mod(surface->conic_y, W);
    build_basic_polynomials(surface);
    a_at_1 = mul_mod(mul_mod(R, R), inv3);
    a_at_lambda = mul_mod(mul_mod(surface->q, surface->q), inv3);
    delta_1 = sub_mod(a_at_1, a_at_0);
    delta_lambda = sub_mod(a_at_lambda, a_at_0);
    interpolation_denominator = mul_mod(lambda, sub_mod(lambda, 1));
    surface->a[0] = a_at_0;
    surface->a[2] = mul_mod(
        sub_mod(delta_lambda, mul_mod(lambda, delta_1)),
        inverse_mod(interpolation_denominator));
    surface->a[1] = sub_mod(delta_1, surface->a[2]);
    surface->beta[0] = mul_mod(W, inv3);
    surface->beta[1] = mul_mod(
        mul_mod(sub_mod(mul_mod(R, surface->conic_x), 1), W), inv3);
    surface->gamma[0] = mul_mod(W, W);
    surface->gamma[1] = mul_mod(
        sub_mod(mul_mod(surface->conic_x, surface->conic_x), 1),
        mul_mod(W, W));
    surface->d = mul_mod(surface->a[2], surface->gamma[1]);
    if (surface->gamma[1] == 0 || !relation_holds(surface)) {
        fputs("internal error: split-cover relation failed\n", stderr);
        exit(3);
    }
    if (poly_evaluate(surface->a, POLY_WIDTH, 0) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, 1) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, lambda) == 0) {
        fputs("internal error: split tangent root vanished\n", stderr);
        exit(3);
    }
    build_H(surface, H);
    if (poly_degree(H) != 5
        || polynomial_gcd_degree(surface->D, H) != 0) {
        return 1;
    }
    poly_derivative(H, derivative);
    if (polynomial_gcd_degree(H, derivative) != 0) {
        return 1;
    }
    return 2;
}

static void build_p1_square(const int *u_small, const Surface *surface,
                            int *answer) {
    int u[POLY_WIDTH];
    int u2[POLY_WIDTH];
    int u3[POLY_WIDTH];
    int term[POLY_WIDTH];
    int product[POLY_WIDTH];
    int tm1_squared[POLY_WIDTH];
    poly_copy_width(u_small, P1_U_WIDTH, u);
    poly_multiply(u, u, u2);
    poly_multiply(u2, u, u3);
    poly_multiply(surface->tz, u3, answer);
    poly_multiply(surface->a, u2, term);
    poly_add_scaled(answer, term, 3);
    poly_multiply(surface->tm1, surface->beta, term);
    poly_multiply(term, u, product);
    poly_add_scaled(answer, product, neg_mod(6));
    poly_multiply(surface->tm1, surface->tm1, tm1_squared);
    poly_multiply(tm1_squared, surface->gamma, term);
    poly_add_scaled(answer, term, 1);
}

static void build_p2_square(const int *u_small, const Surface *surface,
                            int *answer) {
    int u[POLY_WIDTH];
    int u2[POLY_WIDTH];
    int u3[POLY_WIDTH];
    int term[POLY_WIDTH];
    int product[POLY_WIDTH];
    poly_copy_width(u_small, P2_U_WIDTH, u);
    poly_multiply(u, u, u2);
    poly_multiply(u2, u, u3);
    poly_multiply(surface->z, u3, answer);
    poly_multiply(surface->a, u2, term);
    poly_add_scaled(answer, term, 3);
    poly_multiply(surface->T, surface->beta, term);
    poly_multiply(term, u, product);
    poly_add_scaled(answer, product, neg_mod(6));
    poly_multiply(surface->T, surface->T, term);
    poly_multiply(term, surface->gamma, product);
    poly_add_scaled(answer, product, 1);
}

static void build_p3_square(const int *x_small, const Surface *surface,
                            int *answer) {
    int x[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int x3[POLY_WIDTH];
    int term[POLY_WIDTH];
    int product[POLY_WIDTH];
    poly_copy_width(x_small, X_WIDTH, x);
    poly_multiply(x, x, x2);
    poly_multiply(x2, x, x3);
    memcpy(answer, x3, sizeof(x3));
    poly_multiply(surface->a, x2, term);
    poly_add_scaled(answer, term, 3);
    poly_multiply(surface->D, surface->beta, term);
    poly_multiply(term, x, product);
    poly_add_scaled(answer, product, neg_mod(6));
    poly_multiply(surface->D2, surface->gamma, term);
    poly_add_scaled(answer, term, 1);
}

static int verify_full_section(const int *x_small, int x_width,
                               const int *y_small, int y_width,
                               const Surface *surface) {
    int x[POLY_WIDTH];
    int y[POLY_WIDTH];
    int y2[POLY_WIDTH];
    int expected[POLY_WIDTH];
    poly_copy_width(x_small, x_width, x);
    poly_copy_width(y_small, y_width, y);
    poly_multiply(y, y, y2);
    build_p3_square(x, surface, expected);
    return memcmp(y2, expected, sizeof(y2)) == 0;
}

static size_t enumerate_p1(const Surface *surface, P1Section *sections,
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
        int u_lambda;
        int target_lambda;
        int centered_zero;
        decode_polynomial(code, P1_U_WIDTH, u);
        ++counts->p1_tests;
        if (poly_evaluate(u, P1_U_WIDTH, 1) == 0) {
            continue;
        }
        centered_zero = add_mod(u[0], surface->s0);
        if (centered_zero == 0) {
            continue;
        }
        build_p1_square(u, surface, square);
        if (!polynomial_square_root(square, 4, root)) {
            continue;
        }
        u_lambda = poly_evaluate(u, P1_U_WIDTH, surface->lambda);
        target_lambda = sub_mod(
            mul_mod(surface->q, u_lambda),
            mul_mod(sub_mod(surface->lambda, 1),
                    surface->s_lambda));
        if (target_lambda == 0) {
            continue;
        }
        if (poly_evaluate(root, P1_V_WIDTH, surface->lambda)
            != target_lambda) {
            negate_polynomial(root, P1_V_WIDTH);
        }
        if (poly_evaluate(root, P1_V_WIDTH, surface->lambda)
            != target_lambda) {
            continue;
        }
        if (root[0] != centered_zero && root[0] != neg_mod(centered_zero)) {
            fputs("internal error: P1 branch at zero failed\n", stderr);
            exit(3);
        }
        multiply_small_by_poly(u, P1_U_WIDTH, surface->tz, full_x);
        multiply_small_by_poly(root, P1_V_WIDTH, surface->tz, full_y);
        if (!verify_full_section(full_x, X_WIDTH, full_y, Y_WIDTH, surface)) {
            fputs("internal error: P1 full section failed\n", stderr);
            exit(3);
        }
        memcpy(sections[found].u, u, sizeof(u));
        memcpy(sections[found].v, root,
               (size_t)P1_V_WIDTH * sizeof(*root));
        audit_tag(101);
        audit_integer(surface->lambda);
        audit_integer(surface->m);
        audit_integer(surface->R);
        audit_integer(surface->W);
        audit_polynomial(u, P1_U_WIDTH);
        audit_polynomial(root, P1_V_WIDTH);
        ++found;
        ++counts->p1_sections;
    }
    return found;
}

static size_t enumerate_p2(const Surface *surface, P2Section *sections,
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
        int u_lambda;
        int target_lambda;
        decode_polynomial(code, P2_U_WIDTH, u);
        ++counts->p2_tests;
        if (u[0] == 0 || poly_evaluate(u, P2_U_WIDTH, 1) == 0) {
            continue;
        }
        build_p2_square(u, surface, square);
        if (!polynomial_square_root(square, 5, root)) {
            continue;
        }
        u_lambda = poly_evaluate(u, P2_U_WIDTH, surface->lambda);
        target_lambda = add_mod(
            neg_mod(mul_mod(surface->q, u_lambda)),
            mul_mod(mul_mod(surface->lambda,
                            sub_mod(surface->lambda, 1)),
                    surface->s_lambda));
        if (target_lambda == 0) {
            continue;
        }
        if (poly_evaluate(root, P2_V_WIDTH, surface->lambda)
            != target_lambda) {
            negate_polynomial(root, P2_V_WIDTH);
        }
        if (poly_evaluate(root, P2_V_WIDTH, surface->lambda)
            != target_lambda) {
            continue;
        }
        multiply_small_by_poly(u, P2_U_WIDTH, surface->z, full_x);
        multiply_small_by_poly(root, P2_V_WIDTH, surface->z, full_y);
        if (!verify_full_section(full_x, X_WIDTH, full_y, Y_WIDTH, surface)) {
            fputs("internal error: P2 full section failed\n", stderr);
            exit(3);
        }
        memcpy(sections[found].u, u, sizeof(u));
        memcpy(sections[found].v, root,
               (size_t)P2_V_WIDTH * sizeof(*root));
        audit_tag(102);
        audit_integer(surface->lambda);
        audit_integer(surface->m);
        audit_integer(surface->R);
        audit_integer(surface->W);
        audit_polynomial(u, P2_U_WIDTH);
        audit_polynomial(root, P2_V_WIDTH);
        ++found;
        ++counts->p2_sections;
    }
    return found;
}

static size_t enumerate_p3(const Surface *surface, P3Section *sections,
                           Counts *counts) {
    uint64_t limit = integer_power_u64(prime, X_WIDTH);
    size_t found = 0U;
    uint64_t code;
    for (code = 0; code < limit; ++code) {
        int x[X_WIDTH];
        int square[POLY_WIDTH];
        int root[POLY_WIDTH];
        int root_is_zero;
        int sign;
        decode_polynomial(code, X_WIDTH, x);
        ++counts->p3_tests;
        if (x[0] == 0 || poly_evaluate(x, X_WIDTH, 1) == 0
            || poly_evaluate(x, X_WIDTH, surface->lambda) == 0) {
            continue;
        }
        build_p3_square(x, surface, square);
        if (!polynomial_square_root(square, 6, root)) {
            continue;
        }
        root_is_zero = poly_degree(root) < 0;
        for (sign = 0; sign < (root_is_zero ? 1 : 2); ++sign) {
            if (sign != 0) {
                negate_polynomial(root, Y_WIDTH);
            }
            if (!verify_full_section(x, X_WIDTH, root, Y_WIDTH, surface)) {
                fputs("internal error: P3 full section failed\n", stderr);
                exit(3);
            }
            memcpy(sections[found].x, x, sizeof(x));
            memcpy(sections[found].y, root,
                   (size_t)Y_WIDTH * sizeof(*root));
            audit_tag(103);
            audit_integer(surface->lambda);
            audit_integer(surface->m);
            audit_integer(surface->R);
            audit_integer(surface->W);
            audit_polynomial(x, X_WIDTH);
            audit_polynomial(root, Y_WIDTH);
            ++found;
            ++counts->p3_sections;
        }
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

static int section_intersection_degree(const int *x1, const int *y1,
                                       const int *x2, const int *y2) {
    int dx[POLY_WIDTH];
    int dy[POLY_WIDTH];
    int i;
    memcpy(dx, x1, sizeof(dx));
    memcpy(dy, y1, sizeof(dy));
    for (i = 0; i < POLY_WIDTH; ++i) {
        dx[i] = sub_mod(dx[i], x2[i]);
        dy[i] = sub_mod(dy[i], y2[i]);
    }
    return polynomial_gcd_degree(dx, dy);
}

static void p1_full(const P1Section *section, const Surface *surface,
                    int *x, int *y) {
    multiply_small_by_poly(section->u, P1_U_WIDTH, surface->tz, x);
    multiply_small_by_poly(section->v, P1_V_WIDTH, surface->tz, y);
}

static void p2_full(const P2Section *section, const Surface *surface,
                    int *x, int *y) {
    multiply_small_by_poly(section->u, P2_U_WIDTH, surface->z, x);
    multiply_small_by_poly(section->v, P2_V_WIDTH, surface->z, y);
}

static int limiting_point_is_smooth(int x4, int y6) {
    return x4 != 0 && y6 != 0
           && mul_mod(y6, y6) == mul_mod(x4, mul_mod(x4, x4));
}

static int limiting_points_differ(int x1, int y1, int x2, int y2) {
    return x1 != x2 || y1 != y2;
}

static int infinity_is_safe(const P1Section *p1, const P2Section *p2,
                            const P3Section *p3,
                            const Surface *surface) {
    int x1[POLY_WIDTH];
    int y1[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int y2[POLY_WIDTH];
    int x3 = p3->x[4];
    int y3 = p3->y[6];
    p1_full(p1, surface, x1, y1);
    p2_full(p2, surface, x2, y2);
    if (!limiting_point_is_smooth(x1[4], y1[6])
        || !limiting_point_is_smooth(x2[4], y2[6])
        || !limiting_point_is_smooth(x3, y3)) {
        return 0;
    }
    return limiting_points_differ(x1[4], y1[6], x2[4], y2[6])
           && limiting_points_differ(x1[4], y1[6], x3, y3)
           && limiting_points_differ(x2[4], y2[6], x3, y3);
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
        fputs("failed to allocate section storage\n", stderr);
        exit(2);
    }
    return answer;
}

static void print_array(const int *values, int width) {
    int i;
    putchar('[');
    for (i = 0; i < width; ++i) {
        printf("%s%d", i == 0 ? "" : ",", values[i]);
    }
    putchar(']');
}

static void print_surface_parameters(const Surface *surface) {
    printf("\"parameters\":{\"lambda\":%d,\"m\":%d,\"R\":%d,"
           "\"W\":%d,\"x\":%d,\"y\":%d,\"q\":%d},"
           "\"a\":",
           surface->lambda, surface->m, surface->R, surface->W,
           surface->conic_x, surface->conic_y, surface->q);
    print_array(surface->a, 3);
    printf(",\"beta\":");
    print_array(surface->beta, 2);
    printf(",\"gamma\":");
    print_array(surface->gamma, 2);
}

static void print_complete_example(const Surface *surface,
                                   const P1Section *p1,
                                   const P2Section *p2,
                                   const P3Section *p3) {
    printf("{\"event\":\"complete_gram_example\",\"prime\":%d,",
           prime);
    print_surface_parameters(surface);
    printf(",\"P1\":{\"U\":");
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

static void audit_surface(const Surface *surface) {
    audit_tag(100);
    audit_integer(surface->lambda);
    audit_integer(surface->m);
    audit_integer(surface->R);
    audit_integer(surface->W);
    audit_integer(surface->conic_x);
    audit_integer(surface->conic_y);
    audit_integer(surface->q);
    audit_polynomial(surface->a, 3);
    audit_polynomial(surface->beta, 2);
    audit_polynomial(surface->gamma, 2);
}

static void audit_triple(int tag, const Surface *surface, size_t i,
                         size_t j, size_t k) {
    audit_tag(tag);
    audit_integer(surface->lambda);
    audit_integer(surface->m);
    audit_integer(surface->R);
    audit_integer(surface->W);
    audit_u64((uint64_t)i);
    audit_u64((uint64_t)j);
    audit_u64((uint64_t)k);
}

int main(int argc, char **argv) {
    Counts counts = {0};
    int max_print = 20;
    int printed = 0;
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
    if (prime > 13) {
        fputs("PRIME above 13 is refused by this exhaustive p^4/p^5 engine\n",
              stderr);
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
    p3_capacity = checked_capacity(
        UINT64_C(2) * integer_power_u64(prime, X_WIDTH),
        sizeof(*p3_sections));
    p1_sections = checked_calloc(p1_capacity, sizeof(*p1_sections));
    p2_sections = checked_calloc(p2_capacity, sizeof(*p2_sections));
    p3_sections = checked_calloc(p3_capacity, sizeof(*p3_sections));
    printf("{\"event\":\"start\",\"prime\":%d,\"max_print\":%d,"
           "\"family\":\"full split II*+3I3 rational chart\","
           "\"cascade\":\"P1 then P2 then eligible P1/P2 pairs then P3\"}\n",
           prime, max_print);
    for (lambda = 2; lambda < prime; ++lambda) {
        int m;
        for (m = 0; m < prime; ++m) {
            int R;
            for (R = 1; R < prime; ++R) {
                int W;
                for (W = 1; W < prime; ++W) {
                    Surface surface;
                    int surface_status;
                    size_t p1_count;
                    size_t p2_count;
                    size_t p3_count;
                    uint64_t eligible_pairs = UINT64_C(0);
                    uint64_t complete_before;
                    size_t i;
                    ++counts.raw_parameter_tuples;
                    surface_status = build_surface(lambda, m, R, W,
                                                   &surface);
                    if (surface_status == 0) {
                        continue;
                    }
                    ++counts.rational_chart_tuples;
                    if (surface_status != 2) {
                        continue;
                    }
                    ++counts.kodaira_surfaces;
                    audit_surface(&surface);
                    p1_count = enumerate_p1(&surface, p1_sections, &counts);
                    if (p1_count == 0U) {
                        continue;
                    }
                    ++counts.p1_stage_surfaces;
                    p2_count = enumerate_p2(&surface, p2_sections, &counts);
                    if (p2_count == 0U) {
                        continue;
                    }
                    ++counts.p2_stage_surfaces;
                    for (i = 0U; i < p1_count; ++i) {
                        size_t j;
                        for (j = 0U; j < p2_count; ++j) {
                            if (p1_p2_intersection_degree(&p1_sections[i],
                                                          &p2_sections[j])
                                == 2) {
                                ++eligible_pairs;
                                ++counts.p1_p2_intersection_two_pairs;
                                audit_tag(150);
                                audit_u64((uint64_t)i);
                                audit_u64((uint64_t)j);
                            }
                        }
                    }
                    if (eligible_pairs == 0U) {
                        continue;
                    }
                    ++counts.p3_stage_surfaces;
                    p3_count = enumerate_p3(&surface, p3_sections, &counts);
                    if (p3_count == 0U) {
                        continue;
                    }
                    ++counts.pre_candidate_surfaces;
                    complete_before = counts.complete_gram_triples;
                    for (i = 0U; i < p1_count; ++i) {
                        size_t j;
                        int x1[POLY_WIDTH];
                        int y1[POLY_WIDTH];
                        p1_full(&p1_sections[i], &surface, x1, y1);
                        for (j = 0U; j < p2_count; ++j) {
                            size_t k;
                            int x2[POLY_WIDTH];
                            int y2[POLY_WIDTH];
                            if (p1_p2_intersection_degree(&p1_sections[i],
                                                          &p2_sections[j])
                                != 2) {
                                continue;
                            }
                            p2_full(&p2_sections[j], &surface, x2, y2);
                            for (k = 0U; k < p3_count; ++k) {
                                int x3[POLY_WIDTH];
                                int y3[POLY_WIDTH];
                                poly_copy_width(p3_sections[k].x, X_WIDTH, x3);
                                poly_copy_width(p3_sections[k].y, Y_WIDTH, y3);
                                if (section_intersection_degree(x1, y1, x3, y3)
                                        != 2
                                    || section_intersection_degree(x2, y2,
                                                                   x3, y3)
                                           != 2) {
                                    continue;
                                }
                                ++counts.affine_gram_triples;
                                audit_triple(201, &surface, i, j, k);
                                if (!infinity_is_safe(&p1_sections[i],
                                                      &p2_sections[j],
                                                      &p3_sections[k],
                                                      &surface)) {
                                    continue;
                                }
                                ++counts.complete_gram_triples;
                                audit_triple(202, &surface, i, j, k);
                                if (printed < max_print) {
                                    print_complete_example(
                                        &surface, &p1_sections[i],
                                        &p2_sections[j], &p3_sections[k]);
                                    ++printed;
                                }
                            }
                        }
                    }
                    if (counts.complete_gram_triples > complete_before) {
                        ++counts.complete_surfaces;
                    }
                }
            }
        }
    }
    printf("{\"event\":\"summary\",\"prime\":%d,"
           "\"raw_parameter_tuples\":%" PRIu64 ","
           "\"rational_chart_tuples\":%" PRIu64 ","
           "\"kodaira_surfaces\":%" PRIu64 ","
           "\"stage_surfaces\":{\"P1\":%" PRIu64
           ",\"P2\":%" PRIu64 ",\"P3\":%" PRIu64 "},"
           "\"tests\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"sections\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"P1_P2_intersection_two_pairs\":%" PRIu64 ","
           "\"pre_candidate_surfaces\":%" PRIu64 ","
           "\"affine_gram_triples\":%" PRIu64 ","
           "\"complete_gram_triples\":%" PRIu64 ","
           "\"complete_surfaces\":%" PRIu64 ","
           "\"printed_complete_examples\":%d,"
           "\"audit_fnv1a64\":\"%016" PRIx64 "\","
           "\"claim_boundary\":\"exhaustive only on the displayed dense "
           "rational split chart and simple local-branch saturation over F_p; "
           "a modular triple would still require characteristic-zero lifting\"}\n",
           prime, counts.raw_parameter_tuples, counts.rational_chart_tuples,
           counts.kodaira_surfaces, counts.p1_stage_surfaces,
           counts.p2_stage_surfaces, counts.p3_stage_surfaces,
           counts.p1_tests, counts.p2_tests, counts.p3_tests,
           counts.p1_sections, counts.p2_sections, counts.p3_sections,
           counts.p1_p2_intersection_two_pairs,
           counts.pre_candidate_surfaces, counts.affine_gram_triples,
           counts.complete_gram_triples, counts.complete_surfaces, printed,
           audit_hash);
    free(p1_sections);
    free(p2_sections);
    free(p3_sections);
    return 0;
}
