/* Deterministic exact random search for a finite-field quadratic-system point.

   Input header: p dimension rows monomials.  Each following row uses the
   monomial order 1, u_i, u_i^2, u_i*u_j (i<j).
*/

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static int evaluate(
    const int *row, const int *values, int dimension, int prime
) {
    int cursor = 0;
    int64_t total = row[cursor++];
    for (int i = 0; i < dimension; ++i) {
        total += (int64_t)row[cursor++] * values[i];
    }
    for (int i = 0; i < dimension; ++i) {
        total += (int64_t)row[cursor++] * values[i] * values[i];
    }
    for (int i = 0; i < dimension; ++i) {
        for (int j = i + 1; j < dimension; ++j) {
            total += (int64_t)row[cursor++] * values[i] * values[j];
        }
    }
    int residue = (int)(total % prime);
    return residue < 0 ? residue + prime : residue;
}

int main(int argc, char **argv) {
    if (argc != 4 && argc != 5) {
        fprintf(stderr, "usage: %s input max_attempts seed [max_solutions]\n", argv[0]);
        return 2;
    }
    FILE *input = fopen(argv[1], "r");
    if (input == NULL) {
        perror("fopen");
        return 2;
    }
    int prime, dimension, rows, monomials;
    if (fscanf(input, "%d %d %d %d", &prime, &dimension, &rows, &monomials) != 4) {
        fprintf(stderr, "invalid header\n");
        return 2;
    }
    int expected = 1 + 2 * dimension + dimension * (dimension - 1) / 2;
    if (monomials != expected || prime <= 2 || dimension <= 0 || rows <= 0) {
        fprintf(stderr, "invalid dimensions\n");
        return 2;
    }
    int *coefficients = malloc((size_t)rows * monomials * sizeof(int));
    int *values = malloc((size_t)dimension * sizeof(int));
    if (coefficients == NULL || values == NULL) {
        fprintf(stderr, "allocation failure\n");
        return 2;
    }
    for (int i = 0; i < rows * monomials; ++i) {
        if (fscanf(input, "%d", &coefficients[i]) != 1) {
            fprintf(stderr, "truncated coefficient input\n");
            return 2;
        }
        coefficients[i] %= prime;
    }
    fclose(input);
    uint64_t max_attempts = strtoull(argv[2], NULL, 10);
    uint64_t state = strtoull(argv[3], NULL, 10);
    uint64_t max_solutions = argc == 5 ? strtoull(argv[4], NULL, 10) : 1;
    uint64_t solutions = 0;
    for (uint64_t attempt = 1; attempt <= max_attempts; ++attempt) {
        for (int i = 0; i < dimension; ++i) {
            values[i] = (int)(splitmix64(&state) % (uint64_t)prime);
        }
        int good = 1;
        for (int row = 0; row < rows; ++row) {
            if (evaluate(
                    coefficients + row * monomials,
                    values,
                    dimension,
                    prime
                ) != 0) {
                good = 0;
                break;
            }
        }
        if (good) {
            printf("{\"status\":\"found\",\"attempt\":%" PRIu64 ",\"point\":[", attempt);
            for (int i = 0; i < dimension; ++i) {
                printf("%s%d", i ? "," : "", values[i]);
            }
            printf("]}\n");
            ++solutions;
            if (solutions >= max_solutions) {
                printf("{\"status\":\"complete\",\"attempts\":%" PRIu64 ",\"solutions\":%" PRIu64 "}\n", attempt, solutions);
                free(coefficients);
                free(values);
                return 0;
            }
        }
    }
    printf("{\"status\":\"complete\",\"attempts\":%" PRIu64 ",\"solutions\":%" PRIu64 "}\n", max_attempts, solutions);
    free(coefficients);
    free(values);
    return solutions ? 0 : 1;
}
