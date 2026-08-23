/*
 * Exhaustive target-Gram search on the three intrinsic boundary charts of
 * the dense two-split II*+3I3 chart:
 *
 *   A. s0_zero(lambda,q,S),
 *   B. s_lambda_zero(lambda,q,W),
 *   C. node1_beta_gamma_zero(q,W,h), lambda=1-q^2.
 *
 * Exact formulas are certified independently by
 * derive_e8_a2_semistable_two_split_boundaries.py.  Together with the dense
 * (lambda,y,z,W) chart, these cover the marked two-split target locus before
 * the residual-quintic Kodaira open is imposed.
 *
 * This program reuses the audited polynomial and section engine from
 * search_e8_a2_semistable_target.c.  Surface.s0 and Surface.s_lambda carry
 * the exact local node values, including zero on the relevant boundaries,
 * so the simple component-branch tests remain valid without division.
 *
 * A complete triple certifies the target Gram over F_p only.  Absence of a
 * triple is a bounded finite-field result, not nonexistence over Q.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wconversion -Wshadow -pedantic -Werror \
 *      -o /tmp/search_e8_a2_semistable_two_split_boundaries \
 *      research/search_e8_a2_semistable_two_split_boundaries.c
 *
 * Run:
 *   /tmp/search_e8_a2_semistable_two_split_boundaries 5
 *   /tmp/search_e8_a2_semistable_two_split_boundaries 7 20
 */

#define build_surface build_dense_full_split_surface
#define main dense_full_split_engine_main
#define print_complete_example print_dense_full_split_complete_example
#define print_surface_parameters print_dense_full_split_surface_parameters
#include "search_e8_a2_semistable_target.c"
#undef print_surface_parameters
#undef print_complete_example
#undef main
#undef build_surface

typedef enum {
    CHART_S0_ZERO = 0,
    CHART_S_LAMBDA_ZERO = 1,
    CHART_NODE1_ZERO = 2
} ChartKind;

typedef struct {
    ChartKind kind;
    const char *name;
    Counts counts;
    uint64_t pair_gcd_degrees[8];
    int max_print;
    int joint_printed;
    int complete_printed;
} BoundaryContext;

static void interpolate_a(Surface *surface, int value_0, int value_1,
                          int value_lambda) {
    int delta_1 = sub_mod(value_1, value_0);
    int delta_lambda = sub_mod(value_lambda, value_0);
    int denominator = mul_mod(surface->lambda,
                              sub_mod(surface->lambda, 1));
    surface->a[0] = value_0;
    surface->a[2] = mul_mod(
        sub_mod(delta_lambda, mul_mod(surface->lambda, delta_1)),
        inverse_mod(denominator));
    surface->a[1] = sub_mod(delta_1, surface->a[2]);
}

static int finish_boundary_surface(Surface *surface) {
    int H[POLY_WIDTH];
    int derivative[POLY_WIDTH];
    if (!relation_holds(surface)) {
        fputs("internal error: boundary relation failed\n", stderr);
        exit(3);
    }
    if (surface->gamma[1] == 0) {
        fputs("internal error: boundary gamma1 vanished\n", stderr);
        exit(3);
    }
    if (poly_evaluate(surface->a, POLY_WIDTH, 0) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, 1) == 0
        || poly_evaluate(surface->a, POLY_WIDTH, surface->lambda) == 0) {
        fputs("internal error: boundary nodal coefficient vanished\n",
              stderr);
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

static int build_s0_zero(int lambda, int q, int S, Surface *surface) {
    int inv3 = inverse_mod(3);
    int inv_lambda = inverse_mod(lambda);
    memset(surface, 0, sizeof(*surface));
    surface->lambda = lambda;
    surface->m = q;
    surface->R = S;
    surface->W = S;
    surface->q = q;
    surface->s0 = 0;
    surface->s_lambda = S;
    surface->conic_x = 0;
    surface->conic_y = S;
    build_basic_polynomials(surface);
    interpolate_a(surface, inv3,
                  mul_mod(mul_mod(mul_mod(q, q), inv3), inv_lambda),
                  mul_mod(mul_mod(q, q), inv3));
    surface->beta[0] = 0;
    surface->beta[1] = mul_mod(
        mul_mod(mul_mod(q, S), inv3), inv_lambda);
    surface->gamma[0] = 0;
    surface->gamma[1] = mul_mod(mul_mod(S, S), inv_lambda);
    surface->d = mul_mod(
        mul_mod(mul_mod(S, S), inv3),
        mul_mod(inv_lambda, inv_lambda));
    return finish_boundary_surface(surface);
}

static int build_s_lambda_zero(int lambda, int q, int W,
                               Surface *surface) {
    int inv3 = inverse_mod(3);
    int inv_lambda = inverse_mod(lambda);
    int lambda_minus_one = sub_mod(lambda, 1);
    int W2 = mul_mod(W, W);
    memset(surface, 0, sizeof(*surface));
    surface->lambda = lambda;
    surface->m = q;
    surface->R = W;
    surface->W = W;
    surface->q = q;
    surface->s0 = W;
    surface->s_lambda = 0;
    surface->conic_x = W;
    surface->conic_y = 0;
    build_basic_polynomials(surface);
    interpolate_a(
        surface, inv3,
        mul_mod(mul_mod(lambda_minus_one, inv3), inv_lambda),
        mul_mod(mul_mod(q, q), inv3));
    surface->beta[0] = mul_mod(W, inv3);
    surface->beta[1] = neg_mod(mul_mod(mul_mod(W, inv3), inv_lambda));
    surface->gamma[0] = W2;
    surface->gamma[1] = neg_mod(mul_mod(W2, inv_lambda));
    surface->d = neg_mod(mul_mod(
        mul_mod(mul_mod(W2, mul_mod(q, q)), inv3),
        mul_mod(mul_mod(inv_lambda, inv_lambda),
                inverse_mod(lambda_minus_one))));
    return finish_boundary_surface(surface);
}

static int build_node1_zero(int q, int W, int h, Surface *surface) {
    int inv3 = inverse_mod(3);
    int q2 = mul_mod(q, q);
    int lambda = sub_mod(1, q2);
    int W2 = mul_mod(W, W);
    memset(surface, 0, sizeof(*surface));
    surface->lambda = lambda;
    surface->m = q;
    surface->R = h;
    surface->W = W;
    surface->q = q;
    surface->s0 = W;
    surface->s_lambda = mul_mod(q, W);
    surface->conic_x = W;
    surface->conic_y = surface->s_lambda;
    build_basic_polynomials(surface);
    interpolate_a(surface, inv3, mul_mod(h, inv3), mul_mod(q2, inv3));
    surface->beta[0] = mul_mod(W, inv3);
    surface->beta[1] = neg_mod(mul_mod(W, inv3));
    surface->gamma[0] = W2;
    surface->gamma[1] = neg_mod(W2);
    surface->d = neg_mod(mul_mod(
        mul_mod(mul_mod(W2, h), inv3), inverse_mod(q2)));
    return finish_boundary_surface(surface);
}

static void print_boundary_parameters(const BoundaryContext *context,
                                      const Surface *surface) {
    if (context->kind == CHART_S0_ZERO) {
        printf("\"chart\":\"s0_zero\",\"parameters\":{\"lambda\":%d,"
               "\"q\":%d,\"S\":%d},",
               surface->lambda, surface->q, surface->R);
    } else if (context->kind == CHART_S_LAMBDA_ZERO) {
        printf("\"chart\":\"s_lambda_zero\",\"parameters\":{"
               "\"lambda\":%d,\"q\":%d,\"W\":%d},",
               surface->lambda, surface->q, surface->W);
    } else {
        printf("\"chart\":\"node1_beta_gamma_zero\",\"parameters\":{"
               "\"q\":%d,\"W\":%d,\"h\":%d,\"lambda\":%d},",
               surface->q, surface->W, surface->R, surface->lambda);
    }
    printf("\"a\":");
    print_array(surface->a, 3);
    printf(",\"beta\":");
    print_array(surface->beta, 2);
    printf(",\"gamma\":");
    print_array(surface->gamma, 2);
}

static void print_boundary_joint(const BoundaryContext *context,
                                 const Surface *surface,
                                 const P1Section *p1,
                                 const P2Section *p2,
                                 int intersection_degree) {
    printf("{\"event\":\"boundary_joint_profile_example\",\"prime\":%d,",
           prime);
    print_boundary_parameters(context, surface);
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

static void print_boundary_complete(const BoundaryContext *context,
                                    const Surface *surface,
                                    const P1Section *p1,
                                    const P2Section *p2,
                                    const P3Section *p3) {
    printf("{\"event\":\"boundary_complete_gram_example\",\"prime\":%d,",
           prime);
    print_boundary_parameters(context, surface);
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

static void process_boundary_surface(BoundaryContext *context,
                                     Surface *surface,
                                     P1Section *p1_sections,
                                     P2Section *p2_sections,
                                     P3Section *p3_sections) {
    size_t p1_count;
    size_t p2_count;
    size_t p3_count;
    uint64_t eligible_pairs = UINT64_C(0);
    uint64_t complete_before;
    size_t i;
    ++context->counts.kodaira_surfaces;
    audit_tag(300 + (int)context->kind);
    audit_surface(surface);
    p1_count = enumerate_p1(surface, p1_sections, &context->counts);
    if (p1_count == 0U) {
        return;
    }
    ++context->counts.p1_stage_surfaces;
    p2_count = enumerate_p2(surface, p2_sections, &context->counts);
    if (p2_count == 0U) {
        return;
    }
    ++context->counts.p2_stage_surfaces;
    for (i = 0U; i < p1_count; ++i) {
        size_t j;
        for (j = 0U; j < p2_count; ++j) {
            int degree = p1_p2_intersection_degree(&p1_sections[i],
                                                    &p2_sections[j]);
            int bucket = degree + 1;
            if (bucket < 0 || bucket > 6) {
                bucket = 7;
            }
            ++context->pair_gcd_degrees[bucket];
            audit_tag(349);
            audit_u64((uint64_t)i);
            audit_u64((uint64_t)j);
            audit_integer(degree);
            if (context->joint_printed < context->max_print) {
                print_boundary_joint(context, surface, &p1_sections[i],
                                     &p2_sections[j], degree);
                ++context->joint_printed;
            }
            if (degree == 2) {
                ++eligible_pairs;
                ++context->counts.p1_p2_intersection_two_pairs;
                audit_tag(350);
                audit_u64((uint64_t)i);
                audit_u64((uint64_t)j);
            }
        }
    }
    if (eligible_pairs == 0U) {
        return;
    }
    ++context->counts.p3_stage_surfaces;
    p3_count = enumerate_p3(surface, p3_sections, &context->counts);
    if (p3_count == 0U) {
        return;
    }
    ++context->counts.pre_candidate_surfaces;
    complete_before = context->counts.complete_gram_triples;
    for (i = 0U; i < p1_count; ++i) {
        size_t j;
        int x1[POLY_WIDTH];
        int y1[POLY_WIDTH];
        p1_full(&p1_sections[i], surface, x1, y1);
        for (j = 0U; j < p2_count; ++j) {
            size_t k;
            int x2[POLY_WIDTH];
            int y2[POLY_WIDTH];
            if (p1_p2_intersection_degree(&p1_sections[i],
                                          &p2_sections[j]) != 2) {
                continue;
            }
            p2_full(&p2_sections[j], surface, x2, y2);
            for (k = 0U; k < p3_count; ++k) {
                int x3[POLY_WIDTH];
                int y3[POLY_WIDTH];
                poly_copy_width(p3_sections[k].x, X_WIDTH, x3);
                poly_copy_width(p3_sections[k].y, Y_WIDTH, y3);
                if (section_intersection_degree(x1, y1, x3, y3) != 2
                    || section_intersection_degree(x2, y2, x3, y3) != 2) {
                    continue;
                }
                ++context->counts.affine_gram_triples;
                audit_triple(401, surface, i, j, k);
                if (!infinity_is_safe(&p1_sections[i], &p2_sections[j],
                                      &p3_sections[k], surface)) {
                    continue;
                }
                ++context->counts.complete_gram_triples;
                audit_triple(402, surface, i, j, k);
                if (context->complete_printed < context->max_print) {
                    print_boundary_complete(context, surface,
                                            &p1_sections[i],
                                            &p2_sections[j],
                                            &p3_sections[k]);
                    ++context->complete_printed;
                }
            }
        }
    }
    if (context->counts.complete_gram_triples > complete_before) {
        ++context->counts.complete_surfaces;
    }
}

static void print_boundary_summary(const BoundaryContext *context) {
    const Counts *counts = &context->counts;
    printf("{\"event\":\"boundary_chart_summary\",\"prime\":%d,"
           "\"chart\":\"%s\",\"raw_parameter_tuples\":%" PRIu64 ","
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
           "\"audit_fnv1a64\":\"%016" PRIx64 "\","
           "\"claim_boundary\":\"exhaustive only on this displayed "
           "three-dimensional boundary chart and simple local contacts over "
           "F_p; characteristic-zero lifting remains separate\"}\n",
           prime, context->name, counts->raw_parameter_tuples,
           counts->rational_chart_tuples, counts->kodaira_surfaces,
           counts->p1_stage_surfaces, counts->p2_stage_surfaces,
           counts->p3_stage_surfaces, counts->p1_tests, counts->p2_tests,
           counts->p3_tests, counts->p1_sections, counts->p2_sections,
           counts->p3_sections, context->pair_gcd_degrees[0],
           context->pair_gcd_degrees[1], context->pair_gcd_degrees[2],
           context->pair_gcd_degrees[3], context->pair_gcd_degrees[4],
           context->pair_gcd_degrees[5], context->pair_gcd_degrees[6],
           context->pair_gcd_degrees[7],
           counts->p1_p2_intersection_two_pairs,
           counts->pre_candidate_surfaces, counts->affine_gram_triples,
           counts->complete_gram_triples, counts->complete_surfaces,
           audit_hash);
}

static void run_chart(BoundaryContext *context, P1Section *p1_sections,
                      P2Section *p2_sections, P3Section *p3_sections) {
    audit_hash = UINT64_C(1469598103934665603);
    if (context->kind == CHART_S0_ZERO
        || context->kind == CHART_S_LAMBDA_ZERO) {
        int lambda;
        for (lambda = 2; lambda < prime; ++lambda) {
            int q;
            for (q = 1; q < prime; ++q) {
                int value;
                for (value = 1; value < prime; ++value) {
                    Surface surface;
                    int status;
                    ++context->counts.raw_parameter_tuples;
                    if (context->kind == CHART_S0_ZERO) {
                        status = build_s0_zero(lambda, q, value, &surface);
                    } else {
                        status = build_s_lambda_zero(lambda, q, value,
                                                     &surface);
                    }
                    ++context->counts.rational_chart_tuples;
                    if (status != 2) {
                        continue;
                    }
                    process_boundary_surface(context, &surface, p1_sections,
                                             p2_sections, p3_sections);
                }
            }
        }
    } else {
        int q;
        for (q = 1; q < prime; ++q) {
            int W;
            for (W = 1; W < prime; ++W) {
                int h;
                for (h = 1; h < prime; ++h) {
                    Surface surface;
                    int status;
                    ++context->counts.raw_parameter_tuples;
                    if (mul_mod(q, q) == 1) {
                        continue;
                    }
                    status = build_node1_zero(q, W, h, &surface);
                    ++context->counts.rational_chart_tuples;
                    if (status != 2) {
                        continue;
                    }
                    process_boundary_surface(context, &surface, p1_sections,
                                             p2_sections, p3_sections);
                }
            }
        }
    }
    print_boundary_summary(context);
}

int main(int argc, char **argv) {
    int max_print = 20;
    size_t p1_capacity;
    size_t p2_capacity;
    size_t p3_capacity;
    P1Section *p1_sections;
    P2Section *p2_sections;
    P3Section *p3_sections;
    BoundaryContext contexts[3] = {
        {.kind = CHART_S0_ZERO, .name = "s0_zero"},
        {.kind = CHART_S_LAMBDA_ZERO, .name = "s_lambda_zero"},
        {.kind = CHART_NODE1_ZERO, .name = "node1_beta_gamma_zero"},
    };
    int index;
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
        fputs("PRIME above 13 is refused by this exhaustive engine\n",
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
           "\"family\":\"three intrinsic two-split boundary charts\","
           "\"charts\":[\"s0_zero\",\"s_lambda_zero\","
           "\"node1_beta_gamma_zero\"]}\n",
           prime, max_print);
    for (index = 0; index < 3; ++index) {
        contexts[index].max_print = max_print;
        run_chart(&contexts[index], p1_sections, p2_sections, p3_sections);
    }
    free(p1_sections);
    free(p2_sections);
    free(p3_sections);
    return 0;
}
