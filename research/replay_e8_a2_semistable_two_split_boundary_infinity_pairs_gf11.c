/*
 * Fail-closed replay of the four GF(11) s_lambda_zero P1/P2 pairs whose
 * target intersection decomposes as one affine intersection plus one at
 * infinity.  The program reconstructs each boundary surface, checks the
 * Kodaira open, verifies both full sections and their local component
 * branches, and exhausts all P3 quartics.
 */

#define build_surface replay_hidden_build_surface
#define main replay_hidden_main
#define print_complete_example replay_hidden_print_complete
#define print_surface_parameters replay_hidden_print_surface
#include "search_e8_a2_semistable_target.c"
#undef print_surface_parameters
#undef print_complete_example
#undef main
#undef build_surface

typedef struct {
    int lambda;
    int q;
    int W;
    int expected_a[3];
    int expected_beta[2];
    int expected_gamma[2];
    P1Section p1;
    P2Section p2;
} ReplayRecord;

static const ReplayRecord replay_records[] = {
    {
        9, 2, 2,
        {4, 3, 10}, {8, 4}, {4, 2},
        {{8, 5, 4}, {1, 3, 8, 10, 3}},
        {{7, 8, 6, 4}, {3, 7, 9, 5, 7, 3}},
    },
    {
        9, 4, 3,
        {4, 6, 7}, {1, 6}, {9, 10},
        {{7, 9, 5}, {1, 7, 6, 8, 9}},
        {{1, 2, 2, 5}, {6, 8, 9, 0, 10, 9}},
    },
    {
        9, 7, 3,
        {4, 6, 7}, {1, 6}, {9, 10},
        {{7, 9, 5}, {10, 4, 5, 3, 2}},
        {{1, 2, 2, 5}, {5, 3, 2, 0, 1, 2}},
    },
    {
        9, 9, 2,
        {4, 3, 10}, {8, 4}, {4, 2},
        {{8, 5, 4}, {10, 8, 3, 1, 8}},
        {{7, 8, 6, 4}, {8, 4, 2, 6, 4, 8}},
    },
};

static void replay_interpolate_a(Surface *surface, int value_0,
                                 int value_1, int value_lambda) {
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

static void replay_build_s_lambda_zero(const ReplayRecord *record,
                                       Surface *surface) {
    int inv3 = inverse_mod(3);
    int inv_lambda = inverse_mod(record->lambda);
    int lambda_minus_one = sub_mod(record->lambda, 1);
    int W2 = mul_mod(record->W, record->W);
    int H[POLY_WIDTH];
    int derivative[POLY_WIDTH];
    memset(surface, 0, sizeof(*surface));
    surface->lambda = record->lambda;
    surface->m = record->q;
    surface->R = record->W;
    surface->W = record->W;
    surface->q = record->q;
    surface->s0 = record->W;
    surface->s_lambda = 0;
    build_basic_polynomials(surface);
    replay_interpolate_a(
        surface, inv3,
        mul_mod(mul_mod(lambda_minus_one, inv3), inv_lambda),
        mul_mod(mul_mod(record->q, record->q), inv3));
    surface->beta[0] = mul_mod(record->W, inv3);
    surface->beta[1] = neg_mod(
        mul_mod(mul_mod(record->W, inv3), inv_lambda));
    surface->gamma[0] = W2;
    surface->gamma[1] = neg_mod(mul_mod(W2, inv_lambda));
    surface->d = neg_mod(mul_mod(
        mul_mod(mul_mod(W2, mul_mod(record->q, record->q)), inv3),
        mul_mod(mul_mod(inv_lambda, inv_lambda),
                inverse_mod(lambda_minus_one))));
    if (!relation_holds(surface)) {
        fputs("boundary relation replay failed\n", stderr);
        exit(3);
    }
    if (memcmp(surface->a, record->expected_a,
               sizeof(record->expected_a)) != 0
        || memcmp(surface->beta, record->expected_beta,
                  sizeof(record->expected_beta)) != 0
        || memcmp(surface->gamma, record->expected_gamma,
                  sizeof(record->expected_gamma)) != 0) {
        fputs("boundary surface reconstruction mismatch\n", stderr);
        exit(3);
    }
    build_H(surface, H);
    poly_derivative(H, derivative);
    if (surface->gamma[1] == 0 || poly_degree(H) != 5
        || polynomial_gcd_degree(surface->D, H) != 0
        || polynomial_gcd_degree(H, derivative) != 0) {
        fputs("boundary surface left the Kodaira open\n", stderr);
        exit(3);
    }
}

static int replay_reversed_order(const int *left, const int *right,
                                 int top_degree) {
    int coefficient;
    for (coefficient = top_degree; coefficient >= 0; --coefficient) {
        if (left[coefficient] != right[coefficient]) {
            return top_degree - coefficient;
        }
    }
    return POLY_WIDTH;
}

static int replay_infinity_degree(const int *x1, const int *y1,
                                  const int *x2, const int *y2) {
    int x_order = replay_reversed_order(x1, x2, 4);
    int y_order = replay_reversed_order(y1, y2, 6);
    return x_order < y_order ? x_order : y_order;
}

static int replay_total_degree(const int *x1, const int *y1,
                               const int *x2, const int *y2) {
    return section_intersection_degree(x1, y1, x2, y2)
           + replay_infinity_degree(x1, y1, x2, y2);
}

static void replay_verify_local_profiles(const ReplayRecord *record,
                                         const Surface *surface) {
    int p1_u_lambda = poly_evaluate(record->p1.u, P1_U_WIDTH,
                                    surface->lambda);
    int p2_u_lambda = poly_evaluate(record->p2.u, P2_U_WIDTH,
                                    surface->lambda);
    int p1_target_lambda = mul_mod(surface->q, p1_u_lambda);
    int p2_target_lambda = neg_mod(mul_mod(surface->q, p2_u_lambda));
    if (poly_evaluate(record->p1.u, P1_U_WIDTH, 1) == 0
        || add_mod(record->p1.u[0], surface->s0) == 0
        || p1_target_lambda == 0
        || poly_evaluate(record->p1.v, P1_V_WIDTH, surface->lambda)
               != p1_target_lambda
        || record->p2.u[0] == 0
        || poly_evaluate(record->p2.u, P2_U_WIDTH, 1) == 0
        || p2_target_lambda == 0
        || poly_evaluate(record->p2.v, P2_V_WIDTH, surface->lambda)
               != p2_target_lambda) {
        fputs("local component-profile replay failed\n", stderr);
        exit(3);
    }
}

int main(void) {
    Counts counts = {0};
    size_t p3_capacity;
    P3Section *p3_sections;
    uint64_t total_gram_triples = UINT64_C(0);
    size_t record_index;
    prime = 11;
    p3_capacity = checked_capacity(
        UINT64_C(2) * integer_power_u64(prime, X_WIDTH),
        sizeof(*p3_sections));
    p3_sections = checked_calloc(p3_capacity, sizeof(*p3_sections));
    for (record_index = 0U;
         record_index < sizeof(replay_records) / sizeof(replay_records[0]);
         ++record_index) {
        const ReplayRecord *record = &replay_records[record_index];
        Surface surface;
        int x1[POLY_WIDTH];
        int y1[POLY_WIDTH];
        int x2[POLY_WIDTH];
        int y2[POLY_WIDTH];
        int affine_degree;
        int infinity_degree;
        size_t p3_count;
        size_t k;
        replay_build_s_lambda_zero(record, &surface);
        replay_verify_local_profiles(record, &surface);
        p1_full(&record->p1, &surface, x1, y1);
        p2_full(&record->p2, &surface, x2, y2);
        if (!verify_full_section(x1, X_WIDTH, y1, Y_WIDTH, &surface)
            || !verify_full_section(x2, X_WIDTH, y2, Y_WIDTH, &surface)) {
            fputs("full P1/P2 section replay failed\n", stderr);
            return 3;
        }
        affine_degree = p1_p2_intersection_degree(
            &record->p1, &record->p2);
        infinity_degree = replay_infinity_degree(x1, y1, x2, y2);
        if (affine_degree != 1 || infinity_degree != 1) {
            fputs("boundary pair is not on the (1,1) stratum\n", stderr);
            return 3;
        }
        p3_count = enumerate_p3(&surface, p3_sections, &counts);
        for (k = 0U; k < p3_count; ++k) {
            int x3[POLY_WIDTH];
            int y3[POLY_WIDTH];
            poly_copy_width(p3_sections[k].x, X_WIDTH, x3);
            poly_copy_width(p3_sections[k].y, Y_WIDTH, y3);
            if (replay_total_degree(x1, y1, x3, y3) == 2
                && replay_total_degree(x2, y2, x3, y3) == 2) {
                ++total_gram_triples;
            }
        }
        printf("{\"event\":\"verified_boundary_infinity_pair\","
               "\"index\":%zu,\"lambda\":%d,\"q\":%d,\"W\":%d,"
               "\"affine_degree\":1,\"infinity_degree\":1,"
               "\"p3_sections\":%zu}\n",
               record_index, record->lambda, record->q, record->W,
               p3_count);
    }
    printf("{\"event\":\"summary\",\"prime\":11,"
           "\"records_verified\":4,\"new_collision_pairs\":4,"
           "\"p3_tests\":%" PRIu64 ",\"p3_sections\":%" PRIu64
           ",\"total_gram_triples\":%" PRIu64 "}\n",
           counts.p3_tests, counts.p3_sections, total_gram_triples);
    free(p3_sections);
    return 0;
}
