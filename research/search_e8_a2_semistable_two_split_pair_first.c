/*
 * Deterministic pair-first search beyond the small-prime exhaustive range.
 *
 * The search domain is the dense two-split chart (lambda,y,z,W).  Rather than
 * exhausting every U1 on every surface, an affine permutation samples without
 * replacement from the finite incidence universe
 *
 *   {(lambda,y,z,W,U1) : lambda!=0,1; y,z,W!=0; deg U1<=2}.
 *
 * Every retained P1 triggers an exhaustive P2 search on that surface.  Every
 * P1/P2 pair with total resolved intersection two triggers an exhaustive,
 * streaming P3 search (no p^5 storage).  Here "total" is the affine gcd degree
 * plus the local intersection order at infinity.  Thus P2 and P3 conclusions
 * are exact conditional on the sampled P1 incidence; only the first stage is
 * sampled unless BUDGET equals the full incidence-universe size.
 *
 * Sampling is reproducible: SEED determines an affine permutation modulo the
 * universe size, and the first BUDGET entries are visited.  A deterministic
 * FNV-1a audit covers every sampled index, chart status, retained section and
 * tested target pair.  Any complete hit is printed in full.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wconversion -Wshadow -pedantic -Werror \
 *      -o /tmp/search_e8_a2_semistable_two_split_pair_first \
 *      research/search_e8_a2_semistable_two_split_pair_first.c
 *
 * Run:
 *   /tmp/search_e8_a2_semistable_two_split_pair_first PRIME BUDGET SEED [MAX_PRINT]
 */

#define build_surface build_small_prime_full_split_surface
#define main small_prime_full_split_engine_main
#define print_complete_example print_small_prime_complete_example
#define print_surface_parameters print_small_prime_surface_parameters
#include "search_e8_a2_semistable_target.c"
#undef print_surface_parameters
#undef print_complete_example
#undef main
#undef build_surface

typedef struct {
    uint64_t sampled_incidences;
    uint64_t rational_chart_incidences;
    uint64_t kodaira_incidences;
    uint64_t p1_hits;
    uint64_t p2_nonempty_incidences;
    uint64_t eligible_pairs;
    uint64_t affine_degree_two_pairs;
    uint64_t total_degree_two_with_infinity;
    uint64_t pair_gcd_degrees[8];
    uint64_t pair_infinity_degrees[8];
    uint64_t pair_total_degrees[8];
    Counts section_counts;
    int max_print;
    int pair_printed;
    int complete_printed;
} PairFirstCounts;

static int degree_bucket(int degree) {
    int bucket = degree + 1;
    if (bucket < 0 || bucket > 6) {
        bucket = 7;
    }
    return bucket;
}

static void print_degree_counts(const uint64_t *counts) {
    printf("{\"minus_one\":%" PRIu64 ",\"0\":%" PRIu64
           ",\"1\":%" PRIu64 ",\"2\":%" PRIu64
           ",\"3\":%" PRIu64 ",\"4\":%" PRIu64
           ",\"5\":%" PRIu64 ",\"other\":%" PRIu64 "}",
           counts[0], counts[1], counts[2], counts[3], counts[4],
           counts[5], counts[6], counts[7]);
}

static uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b != UINT64_C(0)) {
        uint64_t remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t value;
    *state += UINT64_C(0x9e3779b97f4a7c15);
    value = *state;
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static int build_pair_first_surface(int lambda, int y, int z_ratio, int W,
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
        fputs("internal error: pair-first chart relation failed\n", stderr);
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

static int test_one_p1(const Surface *surface, const int *u,
                       P1Section *section, Counts *counts) {
    int square[POLY_WIDTH];
    int root[POLY_WIDTH];
    int full_x[POLY_WIDTH];
    int full_y[POLY_WIDTH];
    int u_lambda;
    int target_lambda;
    int centered_zero;
    ++counts->p1_tests;
    if (poly_evaluate(u, P1_U_WIDTH, 1) == 0) {
        return 0;
    }
    centered_zero = add_mod(u[0], surface->s0);
    if (centered_zero == 0) {
        return 0;
    }
    build_p1_square(u, surface, square);
    if (!polynomial_square_root(square, 4, root)) {
        return 0;
    }
    u_lambda = poly_evaluate(u, P1_U_WIDTH, surface->lambda);
    target_lambda = sub_mod(
        mul_mod(surface->q, u_lambda),
        mul_mod(sub_mod(surface->lambda, 1), surface->s_lambda));
    if (target_lambda == 0) {
        return 0;
    }
    if (poly_evaluate(root, P1_V_WIDTH, surface->lambda)
        != target_lambda) {
        negate_polynomial(root, P1_V_WIDTH);
    }
    if (poly_evaluate(root, P1_V_WIDTH, surface->lambda)
            != target_lambda
        || (root[0] != centered_zero
            && root[0] != neg_mod(centered_zero))) {
        return 0;
    }
    multiply_small_by_poly(u, P1_U_WIDTH, surface->tz, full_x);
    multiply_small_by_poly(root, P1_V_WIDTH, surface->tz, full_y);
    if (!verify_full_section(full_x, X_WIDTH, full_y, Y_WIDTH, surface)) {
        fputs("internal error: sampled P1 full identity failed\n", stderr);
        exit(3);
    }
    memcpy(section->u, u, (size_t)P1_U_WIDTH * sizeof(*u));
    memcpy(section->v, root, (size_t)P1_V_WIDTH * sizeof(*root));
    ++counts->p1_sections;
    audit_tag(601);
    audit_integer(surface->lambda);
    audit_integer(surface->m);
    audit_integer(surface->R);
    audit_integer(surface->W);
    audit_polynomial(u, P1_U_WIDTH);
    audit_polynomial(root, P1_V_WIDTH);
    return 1;
}

static void print_pair_first_parameters(const Surface *surface) {
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
                                   const int *x2, const int *y2,
                                   int *affine, int *at_infinity) {
    *affine = section_intersection_degree(x1, y1, x2, y2);
    *at_infinity = infinity_intersection_degree(x1, y1, x2, y2);
    return *affine + *at_infinity;
}

static void print_pair_first_seed(const Surface *surface,
                                  const P1Section *p1,
                                  const P2Section *p2, int affine,
                                  int at_infinity, int total) {
    printf("{\"event\":\"pair_first_target_seed\",\"prime\":%d,",
           prime);
    print_pair_first_parameters(surface);
    printf(",\"P1\":{\"U\":");
    print_array(p1->u, P1_U_WIDTH);
    printf(",\"V\":");
    print_array(p1->v, P1_V_WIDTH);
    printf("},\"P2\":{\"U\":");
    print_array(p2->u, P2_U_WIDTH);
    printf(",\"V\":");
    print_array(p2->v, P2_V_WIDTH);
    printf("},\"P1_P2_intersection\":{\"affine\":%d,"
           "\"at_infinity\":%d,\"total\":%d}}\n",
           affine, at_infinity, total);
}

static void print_pair_first_complete(const Surface *surface,
                                      const P1Section *p1,
                                      const P2Section *p2,
                                      const P3Section *p3,
                                      const int *affine,
                                      const int *at_infinity) {
    printf("{\"event\":\"pair_first_complete_gram_hit\",\"prime\":%d,",
           prime);
    print_pair_first_parameters(surface);
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
    printf("},\"affine_intersections\":[%d,%d,%d],"
           "\"infinity_intersections\":[%d,%d,%d],"
           "\"total_intersections\":[2,2,2],"
           "\"gram_numerator_over_3\":[[8,-1,0],[-1,10,0],[0,0,12]]}\n",
           affine[0], affine[1], affine[2], at_infinity[0],
           at_infinity[1], at_infinity[2]);
}

static void stream_p3_for_pair(const Surface *surface, const P1Section *p1,
                               const P2Section *p2,
                               int p1_p2_affine, int p1_p2_infinity,
                               PairFirstCounts *search_counts) {
    uint64_t limit = integer_power_u64(prime, X_WIDTH);
    uint64_t code;
    int x1[POLY_WIDTH];
    int y1[POLY_WIDTH];
    int x2[POLY_WIDTH];
    int y2[POLY_WIDTH];
    p1_full(p1, surface, x1, y1);
    p2_full(p2, surface, x2, y2);
    for (code = UINT64_C(0); code < limit; ++code) {
        int x[X_WIDTH];
        int square[POLY_WIDTH];
        int root[POLY_WIDTH];
        int root_is_zero;
        int sign;
        decode_polynomial(code, X_WIDTH, x);
        ++search_counts->section_counts.p3_tests;
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
            P3Section p3;
            int x3[POLY_WIDTH];
            int y3[POLY_WIDTH];
            int affine[3];
            int at_infinity[3];
            if (sign != 0) {
                negate_polynomial(root, Y_WIDTH);
            }
            memcpy(p3.x, x, sizeof(p3.x));
            memcpy(p3.y, root, sizeof(p3.y));
            ++search_counts->section_counts.p3_sections;
            audit_tag(603);
            audit_polynomial(x, X_WIDTH);
            audit_polynomial(root, Y_WIDTH);
            poly_copy_width(x, X_WIDTH, x3);
            poly_copy_width(root, Y_WIDTH, y3);
            affine[0] = p1_p2_affine;
            at_infinity[0] = p1_p2_infinity;
            if (full_total_intersection(x1, y1, x3, y3, &affine[1],
                                        &at_infinity[1]) != 2
                || full_total_intersection(x2, y2, x3, y3, &affine[2],
                                           &at_infinity[2]) != 2) {
                continue;
            }
            if (affine[0] == 2 && affine[1] == 2 && affine[2] == 2) {
                ++search_counts->section_counts.affine_gram_triples;
            }
            ++search_counts->section_counts.complete_gram_triples;
            audit_tag(702);
            audit_u64(code);
            audit_integer(sign);
            audit_integer(affine[1]);
            audit_integer(at_infinity[1]);
            audit_integer(affine[2]);
            audit_integer(at_infinity[2]);
            if (search_counts->complete_printed < search_counts->max_print) {
                print_pair_first_complete(surface, p1, p2, &p3, affine,
                                          at_infinity);
                ++search_counts->complete_printed;
            }
        }
    }
}

static void decode_incidence(uint64_t incidence, uint64_t u1_count,
                             int *lambda, int *y, int *z_ratio, int *W,
                             int *u) {
    uint64_t surface_code = incidence / u1_count;
    uint64_t u_code = incidence % u1_count;
    uint64_t nonzero_count = (uint64_t)(prime - 1);
    decode_polynomial(u_code, P1_U_WIDTH, u);
    *W = 1 + (int)(surface_code % nonzero_count);
    surface_code /= nonzero_count;
    *z_ratio = 1 + (int)(surface_code % nonzero_count);
    surface_code /= nonzero_count;
    *y = 1 + (int)(surface_code % nonzero_count);
    surface_code /= nonzero_count;
    *lambda = 2 + (int)surface_code;
}

int main(int argc, char **argv) {
    PairFirstCounts search_counts = {0};
    uint64_t requested_budget;
    uint64_t budget;
    uint64_t seed;
    uint64_t state;
    uint64_t raw_surfaces;
    uint64_t u1_count;
    uint64_t universe;
    uint64_t offset;
    uint64_t step;
    uint64_t incidence;
    uint64_t sample_number;
    size_t p2_capacity;
    P2Section *p2_sections;
    if (argc < 4 || argc > 5) {
        fprintf(stderr, "usage: %s PRIME BUDGET SEED [MAX_PRINT]\n", argv[0]);
        return 2;
    }
    prime = atoi(argv[1]);
    if (!is_prime_integer(prime) || prime <= 3 || prime > 23) {
        fputs("PRIME must be a prime in [5,23]\n", stderr);
        return 2;
    }
    requested_budget = strtoull(argv[2], NULL, 10);
    seed = strtoull(argv[3], NULL, 10);
    if (requested_budget == UINT64_C(0)) {
        fputs("BUDGET must be positive\n", stderr);
        return 2;
    }
    search_counts.max_print = 20;
    if (argc == 5) {
        search_counts.max_print = atoi(argv[4]);
        if (search_counts.max_print < 0) {
            fputs("MAX_PRINT must be nonnegative\n", stderr);
            return 2;
        }
    }
    raw_surfaces = (uint64_t)(prime - 2)
                   * integer_power_u64(prime - 1, 3);
    u1_count = integer_power_u64(prime, P1_U_WIDTH);
    if (raw_surfaces > UINT64_MAX / u1_count) {
        fputs("incidence universe overflows uint64_t\n", stderr);
        return 2;
    }
    universe = raw_surfaces * u1_count;
    budget = requested_budget < universe ? requested_budget : universe;
    state = seed;
    offset = splitmix64(&state) % universe;
    step = splitmix64(&state) % universe;
    if (step == UINT64_C(0)) {
        step = UINT64_C(1);
    }
    while (gcd_u64(step, universe) != UINT64_C(1)) {
        ++step;
        if (step == universe) {
            step = UINT64_C(1);
        }
    }
    p2_capacity = checked_capacity(integer_power_u64(prime, P2_U_WIDTH),
                                   sizeof(*p2_sections));
    p2_sections = checked_calloc(p2_capacity, sizeof(*p2_sections));
    printf("{\"event\":\"start\",\"prime\":%d,"
           "\"requested_budget\":%" PRIu64 ",\"budget\":%" PRIu64
           ",\"seed\":%" PRIu64 ",\"incidence_universe\":%" PRIu64
           ",\"offset\":%" PRIu64 ",\"step\":%" PRIu64
           ",\"sampling\":\"affine permutation without replacement\","
           "\"intersection_mode\":\"affine gcd plus infinity order\","
           "\"exhaustive_first_stage\":%s}\n",
           prime, requested_budget, budget, seed, universe, offset, step,
           budget == universe ? "true" : "false");
    incidence = offset;
    for (sample_number = UINT64_C(0); sample_number < budget;
         ++sample_number) {
        int lambda;
        int y;
        int z_ratio;
        int W;
        int u[P1_U_WIDTH];
        int status;
        Surface surface;
        P1Section p1;
        size_t p2_count;
        size_t j;
        ++search_counts.sampled_incidences;
        decode_incidence(incidence, u1_count, &lambda, &y, &z_ratio, &W, u);
        status = build_pair_first_surface(lambda, y, z_ratio, W, &surface);
        audit_tag(500);
        audit_u64(incidence);
        audit_integer(status);
        if (status != 0) {
            ++search_counts.rational_chart_incidences;
        }
        if (status == 2) {
            ++search_counts.kodaira_incidences;
            if (test_one_p1(&surface, u, &p1,
                            &search_counts.section_counts)) {
                ++search_counts.p1_hits;
                p2_count = enumerate_p2(&surface, p2_sections,
                                        &search_counts.section_counts);
                if (p2_count != 0U) {
                    ++search_counts.p2_nonempty_incidences;
                }
                for (j = 0U; j < p2_count; ++j) {
                    int affine;
                    int at_infinity;
                    int total = p1_p2_total_intersection(
                        &p1, &p2_sections[j], &surface, &affine,
                        &at_infinity);
                    ++search_counts.pair_gcd_degrees[
                        degree_bucket(affine)];
                    ++search_counts.pair_infinity_degrees[
                        degree_bucket(at_infinity)];
                    ++search_counts.pair_total_degrees[
                        degree_bucket(total)];
                    audit_tag(650);
                    audit_u64((uint64_t)j);
                    audit_integer(affine);
                    audit_integer(at_infinity);
                    audit_integer(total);
                    if (affine == 2) {
                        ++search_counts.affine_degree_two_pairs;
                    }
                    if (total != 2) {
                        continue;
                    }
                    ++search_counts.eligible_pairs;
                    if (at_infinity != 0) {
                        ++search_counts.total_degree_two_with_infinity;
                    }
                    if (search_counts.pair_printed < search_counts.max_print) {
                        print_pair_first_seed(&surface, &p1,
                                              &p2_sections[j], affine,
                                              at_infinity, total);
                        ++search_counts.pair_printed;
                    }
                    stream_p3_for_pair(&surface, &p1, &p2_sections[j],
                                       affine, at_infinity,
                                       &search_counts);
                }
            }
        }
        incidence += step;
        if (incidence >= universe) {
            incidence -= universe;
        }
    }
    printf("{\"event\":\"summary\",\"prime\":%d,"
           "\"budget\":%" PRIu64 ",\"seed\":%" PRIu64
           ",\"incidence_universe\":%" PRIu64 ","
           "\"sampled_incidences\":%" PRIu64 ","
           "\"rational_chart_incidences\":%" PRIu64 ","
           "\"kodaira_incidences\":%" PRIu64 ","
           "\"tests\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"sections\":{\"P1\":%" PRIu64 ",\"P2\":%" PRIu64
           ",\"P3\":%" PRIu64 "},"
           "\"P1_hits\":%" PRIu64 ","
           "\"P2_nonempty_incidences\":%" PRIu64 ","
           "\"P1_P2_gcd_degree_counts\":",
           prime, budget, seed, universe, search_counts.sampled_incidences,
           search_counts.rational_chart_incidences,
           search_counts.kodaira_incidences,
           search_counts.section_counts.p1_tests,
           search_counts.section_counts.p2_tests,
           search_counts.section_counts.p3_tests,
           search_counts.section_counts.p1_sections,
           search_counts.section_counts.p2_sections,
           search_counts.section_counts.p3_sections, search_counts.p1_hits,
           search_counts.p2_nonempty_incidences);
    print_degree_counts(search_counts.pair_gcd_degrees);
    printf(",\"P1_P2_infinity_degree_counts\":");
    print_degree_counts(search_counts.pair_infinity_degrees);
    printf(",\"P1_P2_total_degree_counts\":");
    print_degree_counts(search_counts.pair_total_degrees);
    printf(",\"affine_degree_two_pairs\":%" PRIu64
           ",\"eligible_pairs\":%" PRIu64
           ",\"total_degree_two_with_positive_infinity\":%" PRIu64
           ",\"affine_gram_triples\":%" PRIu64
           ",\"complete_gram_triples\":%" PRIu64
           ",\"audit_fnv1a64\":\"%016" PRIx64 "\","
           "\"claim_boundary\":\"P1 incidences are a deterministic "
           "without-replacement sample unless budget equals universe; P2 and "
           "P3 are exhaustive only conditional on each sampled P1 hit; "
           "eligibility and Gram tests use total intersections\"}\n",
           search_counts.affine_degree_two_pairs,
           search_counts.eligible_pairs,
           search_counts.total_degree_two_with_infinity,
           search_counts.section_counts.affine_gram_triples,
           search_counts.section_counts.complete_gram_triples, audit_hash);
    free(p2_sections);
    return 0;
}
