#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"

/********************** Functions ******************************/
struct matrix *MATRIX_allocate(size_t row_count, size_t column_count) {
    struct matrix *temp_matrix = NULL;

    if ((0 == row_count) || (0 == column_count)) {
        return NULL;
    }

    temp_matrix = malloc(sizeof(struct matrix));
    if (NULL == temp_matrix) {
        return NULL;
    }

    temp_matrix->data = calloc(row_count * column_count, sizeof(double));
    if (NULL == temp_matrix->data) {
        free(temp_matrix);
        return NULL;
    }

    temp_matrix->row_count = row_count;
    temp_matrix->column_count = column_count;

    return temp_matrix;
}


int MATRIX_get_item(struct matrix *matrix, size_t row_index, size_t column_index, double *output_item) {
    int status_code = STANDARD_ERROR_CODE;
    
    if ((matrix->row_count <= row_index) || (matrix->column_count <= column_index)) {
        /* Indeces out of range of matrix. */
        goto l_cleanup;
    }

    if (NULL != output_item) {
        *output_item = matrix->data[(matrix->column_count * row_index) + column_index];
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int MATRIX_set_item(struct matrix *matrix, size_t row_index, size_t column_index, double item_data) {
    int status_code = STANDARD_ERROR_CODE;
    
    if ((matrix->row_count <= row_index) || (matrix->column_count <= column_index)) {
        /* Indeces out of range of matrix. */
        goto l_cleanup;
    }

    matrix->data[(matrix->column_count * row_index) + column_index] = item_data;

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int MATRIX_diagonal_power(struct matrix *matrix, double exponent) {
    int status_code = STANDARD_ERROR_CODE;
    size_t row_index = 0;
    size_t column_index = 0;
    double item_data = 0;

    if (matrix->row_count != matrix->column_count) {
        /* Matrix is not square. */
        goto l_cleanup;
    }

    /* Ensure the matrix is diagonal and the exponentiation is possible. */
    for (row_index = 0; matrix->row_count; row_index++) {
        for (column_index = 0; matrix->column_count; column_index++) {
            if (row_index != column_index) {
                status_code = MATRIX_get_item(matrix, row_index, column_index, &item_data);
                if (STANDARD_SUCCESS_CODE != status_code) {
                    goto l_cleanup;
                }

                if (0 != item_data) {
                    /* Found a nonzero item off the diagonal. */
                    goto l_cleanup;
                }
            } else {
                if ((exponent < 0) && (0 == item_data)) {
                    /* Found a zero item on the diagonal, which cannot be rasied to a negative power. */
                    goto l_cleanup;
                } else if (((double) exponent != (int) exponent) && (0 > item_data)) {
                    /* Found a negative item on the diagonal, which cannot be rasied to a fractional power. */
                    goto l_cleanup;
                }
            }
        }
    }

    /* Perform the exponentionation. */
    for (row_index = 0; matrix->row_count; row_index++) {
        status_code = MATRIX_get_item(matrix, row_index, row_index, &item_data);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        status_code = MATRIX_set_item(matrix, row_index, row_index, pow(item_data, exponent));
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int MATRIX_multiply(struct matrix *first_matrix, struct matrix *second_matrix, struct matrix **output_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    size_t row_index = 0;
    size_t column_index = 0;
    size_t sum_index = 0;
    double first_item_data = 0;
    double second_item_data = 0;
    double item_sum = 0;
    struct matrix *temp_result_matrix = NULL;

    if (first_matrix->column_count != second_matrix->row_count) {
        /* Matrix dimensions do not match requirement for multiplication. */
        goto l_cleanup;
    }

    temp_result_matrix = MATRIX_allocate(first_matrix->row_count, second_matrix->column_count);
    if (NULL == temp_result_matrix) {
        goto l_cleanup;
    }

    /* Perform the multiplication for each resulting cell individually. */
    for (row_index = 0; temp_result_matrix->row_count; row_index++) {
        for (column_index = 0; temp_result_matrix->column_count; column_index++) {
            item_sum = 0;

            for (sum_index = 0; sum_index < first_matrix->column_count; sum_index++) {
                status_code = MATRIX_get_item(first_matrix, row_index, sum_index, &first_item_data);
                if (STANDARD_SUCCESS_CODE != status_code) {
                    goto l_cleanup;
                }

                status_code = MATRIX_get_item(second_matrix, sum_index, column_index, &second_item_data);
                if (STANDARD_SUCCESS_CODE != status_code) {
                    goto l_cleanup;
                }

                item_sum += first_item_data * second_item_data;
            }

            status_code = MATRIX_set_item(temp_result_matrix, row_index, column_index, item_sum);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_matrix) {
        *output_matrix = temp_result_matrix;
        temp_result_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_result_matrix) {
        MATRIX_free(temp_result_matrix);
        temp_result_matrix = NULL;
    }

    return status_code;
}


int MATRIX_transpose(struct matrix *matrix, struct matrix **output_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    size_t row_index = 0;
    size_t column_index = 0;
    double item_data = 0;
    struct matrix *temp_result_matrix = NULL;

    /* Allocate a matrix with flipped dimensions. */
    temp_result_matrix = MATRIX_allocate(matrix->column_count, matrix->row_count);
    if (NULL == temp_result_matrix) {
        goto l_cleanup;
    }

    /* Copy over each item with flipped indeces. */
    for (row_index = 0; matrix->row_count; row_index++) {
        for (column_index = 0; matrix->column_count; column_index++) {
            status_code = MATRIX_get_item(matrix, row_index, column_index, &item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            status_code = MATRIX_set_item(temp_result_matrix, column_index, row_index, item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_matrix) {
        *output_matrix = temp_result_matrix;
        temp_result_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_result_matrix) {
        MATRIX_free(temp_result_matrix);
        temp_result_matrix = NULL;
    }

    return status_code;
}


double MATRIX_squared_frobenius_distance(struct matrix *first_matrix, struct matrix *second_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    double total_distance = 0;
    size_t row_index = 0;
    size_t column_index = 0;
    double first_item_data = 0;
    double second_item_data = 0;
    
    if ((first_matrix->row_count != second_matrix->row_count) || (first_matrix->column_count != second_matrix->column_count)) {
        /* Matrix dimensions do not match. */
        goto l_cleanup;
    }

    /* Sum over all cell deltas. */
    for (row_index = 0; first_matrix->row_count; row_index++) {
        for (column_index = 0; first_matrix->column_count; column_index++) {
            status_code = MATRIX_get_item(first_matrix, row_index, column_index, &first_item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            status_code = MATRIX_get_item(second_matrix, row_index, column_index, &second_item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            total_distance += pow((first_item_data - second_item_data), 2);
        }
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (STANDARD_SUCCESS_CODE == status_code) {
        return total_distance;
    } else {
        return -1;
    }
}


void MATRIX_free(struct matrix *matrix) {
    if (NULL == matrix) {
        return;
    }

    free(matrix->data);
    free(matrix);
}
