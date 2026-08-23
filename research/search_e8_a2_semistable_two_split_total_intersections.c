/*
 * Audit the infinity-collision boundary omitted by the simple two-split
 * target engine.  For two polynomial sections P,Q, the exact intersection
 * number is
 *
 *   deg gcd(X_P-X_Q,Y_P-Y_Q) +
 *   min(ord_u(u^4(X_P-X_Q)(1/u)),
 *       ord_u(u^6(Y_P-Y_Q)(1/u)))
 *
 * after the forced opposite-component factor at t=lambda is removed for
 * P1/P2.  This program uses that total, rather than requiring two affine
 * intersections and distinct limiting points.
 */

#define build_surface total_probe_full_split_surface
#define main total_probe_full_split_main
#define print_complete_example total_probe_full_split_complete_example
#define print_surface_parameters total_probe_full_split_surface_parameters
#include "search_e8_a2_semistable_target.c"
#undef print_surface_parameters
#undef print_complete_example
#undef main
#undef build_surface

static int build_two_split_probe_surface(int lambda, int y, int z_ratio,
                                         int W, Surface *surface) {
    int inv3 = inverse_mod(3);
    int one_minus_lambda = sub_mod(1, lambda);
    int g = mul_mod(
        sub_mod(mul_mod(y, y), one_minus_lambda), inverse_mod(lambda));
    int q;
    int a_at_0 = inv3;
    int a_at_1;
    int a_at_lambda;
    int delta_1;
    int delta_lambda;
    int interpolation_denominator;
    int H[POLY_WIDTH];
    int derivative[POLY_WIDTH];
    if (g == 0) {
        return 0;
    }
    q = mul_mod(add_mod(one_minus_lambda, mul_mod(lambda, z_ratio)),
                inverse_mod(y));
    if (q == 0) {
        return 0;
    }
    memset(surface, 0, sizeof(*surface));
    surface->lambda = lambda;
    surface->m = y;
    surface->R = z_ratio;
    surface->W = W;
    surface->conic_x = g;
    surface->conic_y = y;
    surface->q = q;
    surface->s0 = W;
    surface->s_lambda = mul_mod(y, W);
    build_basic_polynomials(surface);
    a_at_1 = mul_mod(mul_mod(mul_mod(z_ratio, z_ratio), inv3),
                     inverse_mod(g));
    a_at_lambda = mul_mod(mul_mod(q, q), inv3);
    delta_1 = sub_mod(a_at_1, a_at_0);
    delta_lambda = sub_mod(a_at_lambda, a_at_0);
    interpolation_denominator = mul_mod(lambda, sub_mod(lambda, 1));
    surface->a[0] = a_at_0;
    surface->a[2] = mul_mod(
        sub_mod(delta_lambda, mul_mod(lambda, delta_1)),
        inverse_mod(interpolation_denominator));
    surface->a[1] = sub_mod(delta_1, surface->a[2]);
    surface->beta[0] = mul_mod(W, inv3);
    surface->beta[1] = mul_mod(mul_mod(sub_mod(z_ratio, 1), W), inv3);
    surface->gamma[0] = mul_mod(W, W);
    surface->gamma[1] = mul_mod(sub_mod(g, 1), mul_mod(W, W));
    surface->d = mul_mod(surface->a[2], surface->gamma[1]);
    if (!relation_holds(surface)) {
        fputs("internal error: two-split relation failed\n", stderr);
        exit(3);
    }
    if (surface->gamma[1] == 0) {
        return 1;
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

static int reversed_difference_order(const int *left, const int *right,
                                     int top_degree) {
    int coefficient;
    for (coefficient = top_degree; coefficient >= 0; --coefficient) {
        if (left[coefficient] != right[coefficient]) {
            return top_degree - coefficient;
        }
    }
    return POLY_WIDTH;
}

static int infinity_intersection_degree(const int *x1, const int *y1,
                                        const int *x2, const int *y2) {
    int x_order = reversed_difference_order(x1, x2, 4);
    int y_order = reversed_difference_order(y1, y2, 6);
    return x_order < y_order ? x_order : y_order;
}

static int p1_p2_total_intersection(const P1Section *p1,
                                    const P2Section *p2,
                                    const Surface *surface,
                                    int *affine, int *at_infinity) {
    int x1[POLY_WIDTH];
    int y1[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int y2[POLY_WIDTH];
    *affine = p1_p2_intersection_degree(p1, p2);
    p1_full(p1, surface, x1, y1);
    p2_full(p2, surface, x2, y2);
    *at_infinity = infinity_intersection_degree(x1, y1, x2, y2);
    return *affine + *at_infinity;
}

static int full_total_intersection(const int *x1, const int *y1,
                                   const int *x2, const int *y2) {
    return section_intersection_degree(x1, y1, x2, y2)
           + infinity_intersection_degree(x1, y1, x2, y2);
}

int main(int argc, char **argv) {
    Counts counts = {0};
    uint64_t pair_strata[4][4] = {{0}};
    uint64_t eligible_pairs = UINT64_C(0);
    uint64_t new_collision_pairs = UINT64_C(0);
    uint64_t p3_tested_surfaces = UINT64_C(0);
    uint64_t total_gram_triples = UINT64_C(0);
    size_t p1_capacity;
    size_t p2_capacity;
    size_t p3_capacity;
    P1Section *p1_sections;
    P2Section *p2_sections;
    P3Section *p3_sections;
    int lambda;
    if (argc != 2) {
        fprintf(stderr, "usage: %s PRIME\n", argv[0]);
        return 2;
    }
    prime = atoi(argv[1]);
    if (!is_prime_integer(prime) || prime <= 3 || prime > 13) {
        fputs("PRIME must be 5, 7, 11 or 13\n", stderr);
        return 2;
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
    for (lambda = 2; lambda < prime; ++lambda) {
        int y;
        for (y = 1; y < prime; ++y) {
            int z_ratio;
            for (z_ratio = 1; z_ratio < prime; ++z_ratio) {
                int W;
                for (W = 1; W < prime; ++W) {
                    Surface surface;
                    size_t p1_count;
                    size_t p2_count;
                    uint64_t pair_count = UINT64_C(0);
                    size_t i;
                    if (build_two_split_probe_surface(
                            lambda, y, z_ratio, W, &surface) != 2) {
                        continue;
                    }
                    p1_count = enumerate_p1(&surface, p1_sections, &counts);
                    if (p1_count == 0U) {
                        continue;
                    }
                    p2_count = enumerate_p2(&surface, p2_sections, &counts);
                    if (p2_count == 0U) {
                        continue;
                    }
                    for (i = 0U; i < p1_count; ++i) {
                        size_t j;
                        for (j = 0U; j < p2_count; ++j) {
                            int affine;
                            int at_infinity;
                            int total = p1_p2_total_intersection(
                                &p1_sections[i], &p2_sections[j], &surface,
                                &affine, &at_infinity);
                            if (affine >= 0 && affine < 4
                                && at_infinity >= 0 && at_infinity < 4) {
                                ++pair_strata[affine][at_infinity];
                            }
                            if (total == 2) {
                                ++pair_count;
                                ++eligible_pairs;
                                if (at_infinity != 0) {
                                    ++new_collision_pairs;
                                }
                            }
                        }
                    }
                    if (pair_count != UINT64_C(0)) {
                        size_t p3_count = enumerate_p3(
                            &surface, p3_sections, &counts);
                        ++p3_tested_surfaces;
                        for (i = 0U; i < p1_count; ++i) {
                            size_t j;
                            for (j = 0U; j < p2_count; ++j) {
                                int affine;
                                int at_infinity;
                                const P1Section *p1 = &p1_sections[i];
                                const P2Section *p2 = &p2_sections[j];
                                int x1[POLY_WIDTH];
                                int y1[POLY_WIDTH];
                                int x2[POLY_WIDTH];
                                int y2[POLY_WIDTH];
                                size_t k;
                                if (p1_p2_total_intersection(
                                        p1, p2, &surface, &affine,
                                        &at_infinity) != 2) {
                                    continue;
                                }
                                p1_full(p1, &surface, x1, y1);
                                p2_full(p2, &surface, x2, y2);
                                for (k = 0U; k < p3_count; ++k) {
                                    int x3[POLY_WIDTH];
                                    int y3[POLY_WIDTH];
                                    poly_copy_width(p3_sections[k].x,
                                                    X_WIDTH, x3);
                                    poly_copy_width(p3_sections[k].y,
                                                    Y_WIDTH, y3);
                                    if (full_total_intersection(
                                            x1, y1, x3, y3) == 2
                                        && full_total_intersection(
                                               x2, y2, x3, y3) == 2) {
                                        ++total_gram_triples;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("{\"prime\":%d,\"pair_strata\":{"
           "\"aff0_inf0\":%" PRIu64 ","
           "\"aff1_inf0\":%" PRIu64 ","
           "\"aff1_inf1\":%" PRIu64 ","
           "\"aff2_inf0\":%" PRIu64 ","
           "\"aff2_inf1\":%" PRIu64 ","
           "\"aff3_inf0\":%" PRIu64 "},"
           "\"eligible_total_two_pairs\":%" PRIu64 ","
           "\"new_collision_pairs\":%" PRIu64 ","
           "\"p3_tested_surfaces\":%" PRIu64 ","
           "\"p3_tests\":%" PRIu64 ","
           "\"p3_sections\":%" PRIu64 ","
           "\"total_gram_triples\":%" PRIu64 "}\n",
           prime, pair_strata[0][0], pair_strata[1][0],
           pair_strata[1][1], pair_strata[2][0], pair_strata[2][1],
           pair_strata[3][0], eligible_pairs, new_collision_pairs,
           p3_tested_surfaces, counts.p3_tests, counts.p3_sections,
           total_gram_triples);
    free(p3_sections);
    free(p2_sections);
    free(p1_sections);
    return 0;
}
