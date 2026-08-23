/*
 * Exhaustive polynomial-section search in the non-isotrivial E8+A2^3 family
 *
 *   E_{k,c,l}: y^2 = x^3 - 3 k^2 f^2 x
 *                      + f^2 (c (t-l)^3 + 2 k^3 f),   f=t(t-1).
 *
 * In characteristic p>3, with c*l*(l-1) nonzero and the residual cubic
 *
 *   c(t-l)^3 + 4 k^3 f
 *
 * squarefree, the geometric fibre configuration is
 *
 *   II* + IV@0 + IV@1 + I3@l + 3 I1.
 *
 * The program exhausts every x in F_p[t] of degree at most four and tests
 * whether the right-hand side is a square y^2 with deg(y)<=6.  It stores one
 * representative from {P,-P}, then classifies the finite A2 component mask:
 *
 *   bit 0: P specializes to the IV cusp (0,0) at t=0;
 *   bit 1: P specializes to the IV cusp (0,0) at t=1;
 *   bit 2: P specializes to the I3 node (k f(l),0) at t=l.
 *
 * A polynomial section with n set bits has Shioda height 4-2n/3.  Thus the
 * diagonal heights 8/3, 10/3, 4 of the target Mordell--Weil Gram require
 * popcounts 2, 1, 0 respectively.  Component masks alone do not certify the
 * off-diagonal pairings or independence.
 *
 * Build:
 *   cc -O3 -Wall -Wextra -pedantic -o /tmp/search_mixed \
 *      research/search_e8_a2_mixed_modular.c
 * Run:
 *   /tmp/search_mixed 7
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH 13
#define XWIDTH 5
#define MAX_STORED_SECTIONS 4096
#define MAX_PRINTED_SURFACES 200

typedef struct {
    int x[XWIDTH];
    int y[7];
    unsigned mask;
} Section;

static int prime;

static int mod_i64(int64_t value) {
    int answer = (int)(value % prime);
    return answer < 0 ? answer + prime : answer;
}

static int add(int a, int b) { return mod_i64((int64_t)a + b); }
static int sub(int a, int b) { return mod_i64((int64_t)a - b); }
static int mul(int a, int b) { return mod_i64((int64_t)a * b); }

static int power_mod(int a, int exponent) {
    int result = 1;
    while (exponent > 0) {
        if (exponent & 1) {
            result = mul(result, a);
        }
        a = mul(a, a);
        exponent >>= 1;
    }
    return result;
}

static int integer_power_checked(int base, int exponent) {
    int64_t result = 1;
    while (exponent-- > 0) {
        result *= base;
        if (result > INT32_MAX) {
            fprintf(stderr, "search space exceeds 32-bit loop counter\n");
            exit(2);
        }
    }
    return (int)result;
}

static int inverse_mod(int a) {
    if (a == 0) {
        fprintf(stderr, "attempted inversion of zero\n");
        exit(2);
    }
    return power_mod(a, prime - 2);
}

static int is_prime_integer(int n) {
    int d;
    if (n < 2) {
        return 0;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (d = 3; (int64_t)d * d <= n; d += 2) {
        if ((n % d) == 0) {
            return 0;
        }
    }
    return 1;
}

static int polynomial_degree(const int *a, int width) {
    int i;
    for (i = width - 1; i >= 0; --i) {
        if (a[i] != 0) {
            return i;
        }
    }
    return -1;
}

static int evaluate(const int *a, int width, int value) {
    int i;
    int result = 0;
    for (i = width - 1; i >= 0; --i) {
        result = add(mul(result, value), a[i]);
    }
    return result;
}

/* Return one square root, if the polynomial is an exact square. */
static int polynomial_square_root(const int *rhs, int y[7]) {
    int degree = polynomial_degree(rhs, WIDTH);
    int half_degree;
    int leading_root = -1;
    int candidate;
    int j;
    int k;

    memset(y, 0, 7 * sizeof(*y));
    if (degree < 0) {
        return 1;
    }
    if (degree & 1) {
        return 0;
    }
    half_degree = degree / 2;
    if (half_degree > 6) {
        return 0;
    }
    for (candidate = 1; candidate < prime; ++candidate) {
        if (mul(candidate, candidate) == rhs[degree]) {
            leading_root = candidate;
            break;
        }
    }
    if (leading_root < 0) {
        return 0;
    }
    y[half_degree] = leading_root;

    for (j = half_degree - 1; j >= 0; --j) {
        int known = 0;
        int denominator = mul(2, y[half_degree]);
        k = half_degree + j;
        for (int i = j + 1; i <= half_degree; ++i) {
            int other = k - i;
            if (other < 0 || other > half_degree || other == j) {
                continue;
            }
            known = add(known, mul(y[i], y[other]));
        }
        y[j] = mul(sub(rhs[k], known), inverse_mod(denominator));
    }

    for (k = 0; k < WIDTH; ++k) {
        int square_coefficient = 0;
        for (int i = 0; i <= 6; ++i) {
            int other = k - i;
            if (other >= 0 && other <= 6) {
                square_coefficient = add(square_coefficient,
                                         mul(y[i], y[other]));
            }
        }
        if (square_coefficient != rhs[k]) {
            return 0;
        }
    }
    return 1;
}

static int cubic_discriminant(int a, int b, int c, int d) {
    /* b^2 c^2 - 4ac^3 - 4b^3d - 27a^2d^2 + 18abcd */
    int result = mul(mul(b, b), mul(c, c));
    result = sub(result, mul(4, mul(a, mul(c, mul(c, c)))));
    result = sub(result, mul(4, mul(mul(b, mul(b, b)), d)));
    result = sub(result, mul(27, mul(mul(a, a), mul(d, d))));
    result = add(result, mul(18, mul(mul(a, b), mul(c, d))));
    return result;
}

static int popcount3(unsigned mask) {
    return (int)(mask & 1U) + (int)((mask >> 1) & 1U)
           + (int)((mask >> 2) & 1U);
}

static void print_section(const Section *section) {
    int i;
    printf("{\"mask\":%u,\"x\":[", section->mask);
    for (i = 0; i < XWIDTH; ++i) {
        printf("%s%d", i == 0 ? "" : ",", section->x[i]);
    }
    printf("],\"y\":[");
    for (i = 0; i < 7; ++i) {
        printf("%s%d", i == 0 ? "" : ",", section->y[i]);
    }
    printf("]}");
}

int main(int argc, char **argv) {
    uint64_t total_x_tests = 0;
    uint64_t total_sections = 0;
    uint64_t generic_surfaces = 0;
    uint64_t nonempty_surfaces = 0;
    uint64_t profile_complete_surfaces = 0;
    int printed_surfaces = 0;

    if (argc != 2) {
        fprintf(stderr, "usage: %s PRIME\n", argv[0]);
        return 2;
    }
    prime = atoi(argv[1]);
    if (!is_prime_integer(prime) || prime <= 3) {
        fprintf(stderr, "PRIME must be a prime integer greater than 3\n");
        return 2;
    }

    printf("{\"event\":\"start\",\"prime\":%d}\n", prime);

    int nonsquare = 2;
    while (nonsquare < prime && power_mod(nonsquare, (prime - 1) / 2) != prime - 1) {
        ++nonsquare;
    }
    if (nonsquare == prime) {
        fprintf(stderr, "failed to locate a quadratic nonsquare\n");
        return 2;
    }

    for (int twist_index = 0; twist_index < 2; ++twist_index) {
      int coefficient_k = twist_index == 0 ? 1 : nonsquare;
      int k2 = mul(coefficient_k, coefficient_k);
      int k3 = mul(k2, coefficient_k);

      for (int lambda = 2; lambda < prime; ++lambda) {
        int lambda2 = mul(lambda, lambda);
        int lambda3 = mul(lambda2, lambda);
        int fibre_f = mul(coefficient_k, mul(lambda, sub(lambda, 1)));

        for (int coefficient_c = 1; coefficient_c < prime; ++coefficient_c) {
            int residual_a = coefficient_c;
            int residual_b = add(mul(-3, mul(coefficient_c, lambda)), mul(4, k3));
            int residual_c = sub(mul(3, mul(coefficient_c, lambda2)), mul(4, k3));
            int residual_d = mul(-1, mul(coefficient_c, lambda3));
            int residual_discriminant = cubic_discriminant(
                residual_a, residual_b, residual_c, residual_d);
            int f2[5] = {0, 0, 1, mod_i64(-2), 1};
            int a_poly[5];
            int inner[4];
            int b_poly[8] = {0};
            Section sections[MAX_STORED_SECTIONS];
            int section_count = 0;
            int profile_counts[4] = {0, 0, 0, 0};

            if (residual_discriminant == 0) {
                continue;
            }
            ++generic_surfaces;

            for (int i = 0; i < 5; ++i) {
                a_poly[i] = mul(mul(-3, k2), f2[i]);
            }
            inner[0] = mul(-1, mul(coefficient_c, lambda3));
            inner[1] = sub(mul(3, mul(coefficient_c, lambda2)), mul(2, k3));
            inner[2] = add(mul(-3, mul(coefficient_c, lambda)), mul(2, k3));
            inner[3] = coefficient_c;
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 4; ++j) {
                    b_poly[i + j] = add(b_poly[i + j],
                                                mul(f2[i], inner[j]));
                }
            }

            for (int encoded_x = 0, limit = integer_power_checked(prime, 5);
                 encoded_x < limit; ++encoded_x) {
                int code = encoded_x;
                int x[XWIDTH];
                int rhs[WIDTH] = {0};
                int y[7];
                unsigned mask = 0;

                for (int i = 0; i < XWIDTH; ++i) {
                    x[i] = code % prime;
                    code /= prime;
                }
                ++total_x_tests;

                for (int i = 0; i < XWIDTH; ++i) {
                    for (int j = 0; j < XWIDTH; ++j) {
                        for (int k = 0; k < XWIDTH; ++k) {
                            rhs[i + j + k] = add(rhs[i + j + k],
                                mul(x[i], mul(x[j], x[k])));
                        }
                    }
                }
                for (int i = 0; i < 5; ++i) {
                    for (int j = 0; j < XWIDTH; ++j) {
                        rhs[i + j] = add(rhs[i + j],
                                        mul(a_poly[i], x[j]));
                    }
                }
                for (int i = 0; i < 8; ++i) {
                    rhs[i] = add(rhs[i], b_poly[i]);
                }

                if (!polynomial_square_root(rhs, y)) {
                    continue;
                }
                if (x[0] == 0 && y[0] == 0) {
                    mask |= 1U;
                }
                if (evaluate(x, XWIDTH, 1) == 0
                    && evaluate(y, 7, 1) == 0) {
                    mask |= 2U;
                }
                if (evaluate(x, XWIDTH, lambda) == fibre_f
                    && evaluate(y, 7, lambda) == 0) {
                    mask |= 4U;
                }
                if (section_count >= MAX_STORED_SECTIONS) {
                    fprintf(stderr, "too many sections on k=%d c=%d lambda=%d\n",
                            coefficient_k, coefficient_c, lambda);
                    return 3;
                }
                memcpy(sections[section_count].x, x, sizeof(x));
                memcpy(sections[section_count].y, y, sizeof(y));
                sections[section_count].mask = mask;
                ++section_count;
                ++profile_counts[popcount3(mask)];
                ++total_sections;
            }

            if (section_count > 0) {
                ++nonempty_surfaces;
                printf("{\"event\":\"sections\",\"prime\":%d,"
                       "\"k\":%d,\"lambda\":%d,\"c\":%d,"
                       "\"residual_discriminant\":%d,"
                       "\"counts_by_nonidentity_fibres\":[%d,%d,%d,%d],"
                       "\"sections\":[",
                       prime, coefficient_k, lambda, coefficient_c,
                       residual_discriminant,
                       profile_counts[0], profile_counts[1],
                       profile_counts[2], profile_counts[3]);
                for (int i = 0; i < section_count; ++i) {
                    printf("%s", i == 0 ? "" : ",");
                    print_section(&sections[i]);
                }
                printf("]}\n");
            }

            if (profile_counts[0] && profile_counts[1] && profile_counts[2]) {
                ++profile_complete_surfaces;
                if (printed_surfaces < MAX_PRINTED_SURFACES) {
                    printf("{\"event\":\"profile_complete\",\"prime\":%d,"
                           "\"k\":%d,\"lambda\":%d,\"c\":%d,"
                           "\"residual_discriminant\":%d,"
                           "\"counts_by_nonidentity_fibres\":[%d,%d,%d,%d],"
                           "\"examples\":[",
                           prime, coefficient_k, lambda, coefficient_c,
                           residual_discriminant,
                           profile_counts[0], profile_counts[1],
                           profile_counts[2], profile_counts[3]);
                    for (int target = 2, emitted = 0; target >= 0; --target) {
                        for (int i = 0; i < section_count; ++i) {
                            if (popcount3(sections[i].mask) == target) {
                                printf("%s", emitted++ == 0 ? "" : ",");
                                print_section(&sections[i]);
                                break;
                            }
                        }
                    }
                    printf("]}\n");
                    ++printed_surfaces;
                }
            }
        }
      }
    }

    printf("{\"event\":\"summary\",\"prime\":%d,"
           "\"generic_surfaces\":%llu,\"nonempty_surfaces\":%llu,"
           "\"x_tests\":%llu,"
           "\"section_pairs\":%llu,"
           "\"profile_complete_surfaces\":%llu,"
           "\"printed_surfaces\":%d}\n",
           prime,
           (unsigned long long)generic_surfaces,
           (unsigned long long)nonempty_surfaces,
           (unsigned long long)total_x_tests,
           (unsigned long long)total_sections,
           (unsigned long long)profile_complete_surfaces,
           printed_surfaces);
    return 0;
}
