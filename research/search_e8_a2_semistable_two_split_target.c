/*
 * Exhaustive target-Gram search on a rational two-split cover of the full
 * II* + 3 I3 semistable E8+A2^3 chart.
 *
 * Only the I3 fibres at 0 and lambda are required to split.  This is the
 * necessary cover for the target profiles: P1 is nonidentity at 0 and
 * lambda, P2 is nonidentity at lambda, while every target section meets the
 * identity component at 1.  Requiring I3@1 to split would discard relevant
 * rational surfaces.
 *
 * Over F_p, p>3, use parameters (lambda,y,z,W).  Put
 *
 *   g = (y^2-(1-lambda))/lambda,
 *   q = ((1-lambda)+lambda*z)/y.
 *
 * Interpolate deg(a,beta,gamma)<=(2,1,1) from
 *
 *   a(0)=1/3,       beta(0)=W/3,       gamma(0)=W^2,
 *   a(1)=z^2/(3g),  beta(1)=zW/3,      gamma(1)=gW^2,
 *   a(lambda)=q^2/3,
 *   beta(lambda)=qyW/3, gamma(lambda)=y^2W^2.
 *
 * The identities qy=(1-lambda)+lambda*z and
 * y^2=(1-lambda)+lambda*g imply
 *
 *   a gamma - 3 beta^2 = dD,  D=t(t-1)(t-lambda).
 *
 * The shared engine below checks this relation coefficient-by-coefficient,
 * checks the exact II*+3I3+5I1 Kodaira open, and exhausts the three section
 * charts, simple local component labels, resolved intersections and safe
 * limiting points at infinity.  A complete triple certifies
 *
 *          1/3 [[8,-1,0],[-1,10,0],[0,0,12]]
 *
 * over F_p only.  A characteristic-zero lift and rank 31 remain separate.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wconversion -Wshadow -pedantic -Werror \
 *      -o /tmp/search_e8_a2_semistable_two_split_target \
 *      research/search_e8_a2_semistable_two_split_target.c
 *
 * Run:
 *   /tmp/search_e8_a2_semistable_two_split_target 5
 *   /tmp/search_e8_a2_semistable_two_split_target 7 20
 *
 * The implementation includes the audited polynomial/section engine from
 * search_e8_a2_semistable_target.c in the same translation unit.  Its chart
 * constructor, formatter and entry point are renamed; the two-split versions
 * below replace them without duplicating the arithmetic core.
 */

#define build_surface build_full_split_surface
#define main full_split_engine_main
#define print_complete_example print_full_split_complete_example
#define print_surface_parameters print_full_split_surface_parameters
#include "search_e8_a2_semistable_target.c"
#undef print_surface_parameters
#undef print_complete_example
#undef main
#undef build_surface

/* Return 0 off this rational chart, 1 on the chart but off the Kodaira open,
 * and 2 on the exact II*+3I3+5I1 open set.  Surface.m stores y, Surface.R
 * stores z, Surface.conic_x stores g, and Surface.conic_y stores y. */
static int build_two_split_surface(int lambda, int y, int z_ratio, int W,
                                   Surface *surface) {
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
    if (poly_evaluate(surface->a, POLY_WIDTH, 0) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, 1) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, lambda) == 0) {
        fputs("internal error: two-split nodal coefficient vanished\n",
              stderr);
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

static void print_two_split_surface_parameters(const Surface *surface) {
    printf("\"parameters\":{\"lambda\":%d,\"y\":%d,\"z\":%d,"
           "\"W\":%d,\"g\":%d,\"q\":%d},\"a\":",
           surface->lambda, surface->m, surface->R, surface->W,
           surface->conic_x, surface->q);
    print_array(surface->a, 3);
    printf(",\"beta\":");
    print_array(surface->beta, 2);
    printf(",\"gamma\":");
    print_array(surface->gamma, 2);
}

static void print_two_split_complete_example(const Surface *surface,
                                             const P1Section *p1,
                                             const P2Section *p2,
                                             const P3Section *p3) {
    printf("{\"event\":\"complete_gram_example\",\"prime\":%d,",
           prime);
    print_two_split_surface_parameters(surface);
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

static void print_two_split_pair_seed(const Surface *surface,
                                      const P1Section *p1,
                                      const P2Section *p2) {
    printf("{\"event\":\"p1_p2_intersection_two_seed\",\"prime\":%d,",
           prime);
    print_two_split_surface_parameters(surface);
    printf(",\"P1\":{\"U\":");
    print_array(p1->u, P1_U_WIDTH);
    printf(",\"V\":");
    print_array(p1->v, P1_V_WIDTH);
    printf("},\"P2\":{\"U\":");
    print_array(p2->u, P2_U_WIDTH);
    printf(",\"V\":");
    print_array(p2->v, P2_V_WIDTH);
    printf("},\"resolved_P1_P2_intersection\":2}\n");
}

static void print_two_split_joint_profile(const Surface *surface,
                                          const P1Section *p1,
                                          const P2Section *p2,
                                          int intersection_degree) {
    printf("{\"event\":\"joint_profile_example\",\"prime\":%d,",
           prime);
    print_two_split_surface_parameters(surface);
    printf(",\"P1\":{\"U\":");
    print_array(p1->u, P1_U_WIDTH);
    printf(",\"V\":");
    print_array(p1->v, P1_V_WIDTH);
    printf("},\"P2\":{\"U\":");
    print_array(p2->u, P2_U_WIDTH);
    printf(",\"V\":");
    print_array(p2->v, P2_V_WIDTH);
    printf("},\"resolved_P1_P2_intersection\":%d,"
           "\"target_intersection_two\":%s}\n",
           intersection_degree,
           intersection_degree == 2 ? "true" : "false");
}

int main(int argc, char **argv) {
    Counts counts = {0};
    int max_print = 20;
    int printed = 0;
    int pair_printed = 0;
    int joint_printed = 0;
    uint64_t pair_gcd_degrees[8] = {0};
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
           "\"family\":\"two-split II*+3I3 rational chart\","
           "\"split_fibres\":[0,\"lambda\"],"
           "\"cascade\":\"P1 then P2 then eligible P1/P2 pairs then P3\"}\n",
           prime, max_print);
    for (lambda = 2; lambda < prime; ++lambda) {
        int y;
        for (y = 1; y < prime; ++y) {
            int z_ratio;
            for (z_ratio = 1; z_ratio < prime; ++z_ratio) {
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
                    surface_status = build_two_split_surface(
                        lambda, y, z_ratio, W, &surface);
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
                            int intersection_degree =
                                p1_p2_intersection_degree(&p1_sections[i],
                                                          &p2_sections[j]);
                            int bucket = intersection_degree + 1;
                            if (bucket < 0 || bucket > 6) {
                                bucket = 7;
                            }
                            ++pair_gcd_degrees[bucket];
                            audit_tag(149);
                            audit_u64((uint64_t)i);
                            audit_u64((uint64_t)j);
                            audit_integer(intersection_degree);
                            if (joint_printed < max_print) {
                                print_two_split_joint_profile(
                                    &surface, &p1_sections[i],
                                    &p2_sections[j], intersection_degree);
                                ++joint_printed;
                            }
                            if (intersection_degree == 2) {
                                ++eligible_pairs;
                                ++counts.p1_p2_intersection_two_pairs;
                                audit_tag(150);
                                audit_u64((uint64_t)i);
                                audit_u64((uint64_t)j);
                                if (pair_printed < max_print) {
                                    print_two_split_pair_seed(
                                        &surface, &p1_sections[i],
                                        &p2_sections[j]);
                                    ++pair_printed;
                                }
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
                                    print_two_split_complete_example(
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
           "\"P1_P2_gcd_degree_counts\":{\"minus_one\":%" PRIu64
           ",\"0\":%" PRIu64 ",\"1\":%" PRIu64
           ",\"2\":%" PRIu64 ",\"3\":%" PRIu64
           ",\"4\":%" PRIu64 ",\"5\":%" PRIu64
           ",\"other\":%" PRIu64 "},"
           "\"P1_P2_intersection_two_pairs\":%" PRIu64 ","
           "\"pre_candidate_surfaces\":%" PRIu64 ","
           "\"affine_gram_triples\":%" PRIu64 ","
           "\"complete_gram_triples\":%" PRIu64 ","
           "\"complete_surfaces\":%" PRIu64 ","
           "\"printed_complete_examples\":%d,"
           "\"audit_fnv1a64\":\"%016" PRIx64 "\","
           "\"claim_boundary\":\"exhaustive only on the displayed dense "
           "two-split chart and simple local-branch saturation over F_p; a "
           "modular triple would still require characteristic-zero lifting\"}\n",
           prime, counts.raw_parameter_tuples, counts.rational_chart_tuples,
           counts.kodaira_surfaces, counts.p1_stage_surfaces,
           counts.p2_stage_surfaces, counts.p3_stage_surfaces,
           counts.p1_tests, counts.p2_tests, counts.p3_tests,
           counts.p1_sections, counts.p2_sections, counts.p3_sections,
           pair_gcd_degrees[0], pair_gcd_degrees[1], pair_gcd_degrees[2],
           pair_gcd_degrees[3], pair_gcd_degrees[4], pair_gcd_degrees[5],
           pair_gcd_degrees[6], pair_gcd_degrees[7],
           counts.p1_p2_intersection_two_pairs,
           counts.pre_candidate_surfaces, counts.affine_gram_triples,
           counts.complete_gram_triples, counts.complete_surfaces, printed,
           audit_hash);
    free(p1_sections);
    free(p2_sections);
    free(p3_sections);
    return 0;
}
