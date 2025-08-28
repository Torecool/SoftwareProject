#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "symnmf.h"
#include "vector.h"
#include "matrix.h"

/********************** Constants *****************************/
#define MAX_ITERATIONS (300)
#define CONVERGENCE_THRESHOLD (1e-4)
#define BETA_VALUE (0.5)

#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)

/********************** SymNMF Algorithm **********************/
/**
 * Calculate the squared Euclidian distance of two given datapoints (of the same dimension).
 * 
 * @param first_point The first datapoint in the calculation.
 * @param second_point The second datapoint in the calculation.
 * @param dimension The dimension of the datapoints.
 * @return The computed distance.
 */
static double symnmf_squared_euclidian_distance(double *first_point, double *second_point, size_t dimension) {
    double temp_distance = 0;
    size_t i = 0;

    for (i = 0; i < dimension; i++) {
        temp_distance += pow(first_point[i] - second_point[i], 2);
    }

    return temp_distance;
}


/**
 * Execute a single iteration of the SymNMF algorithm.
 * In each iteration, we update the H matrix according to a set rule, 
 * using: the previous H matrix, the W matrix and the beta constant.
 * 
 * @param prev_h_matrix The H matrix from the previous iteration.
 * @param w_matrix The W matrix, computed at the start of the algorithm.
 * @param output_next_h_matrix An optional output parameter to contain the new H matrix.
 *  The matrix is dynamically allocated and must be later freed via MATRIX_free().
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int symnmf_update_h_matrix(
    struct matrix *prev_h_matrix,
    struct matrix *w_matrix,
    struct matrix **output_next_h_matrix
) {
    int status_code = STANDARD_ERROR_CODE;
    struct matrix *temp_w_mult_h_matrix = NULL;
    struct matrix *temp_h_transpose_matrix = NULL;
    struct matrix *temp_h_mult_ht_matrix = NULL;
    struct matrix *temp_h_mult_ht_mult_h_matrix = NULL;
    struct matrix *temp_next_h_matrix = NULL;
    size_t row_index = 0;
    size_t column_index = 0;
    double prev_h_item = 0;
    double w_mult_h_item = 0;
    double h_mult_ht_mult_h_item = 0;
    double next_h_item = 0;

    /* Allocate next H(t+1). */
    temp_next_h_matrix = MATRIX_allocate(prev_h_matrix->row_count, prev_h_matrix->column_count);
    if (NULL == temp_next_h_matrix) {
        goto l_cleanup;
    }

    /* Calculate W * H(t). */
    status_code = MATRIX_multiply(w_matrix, prev_h_matrix, &temp_w_mult_h_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Calculate H(t) transpose. */
    status_code = MATRIX_transpose(prev_h_matrix, &temp_h_transpose_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    
    /* Calculate H(t) * H(t)^T. */
    status_code = MATRIX_multiply(prev_h_matrix, temp_h_transpose_matrix, &temp_h_mult_ht_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Calculate (H(t) * H(t)^T) * H(t). */
    status_code = MATRIX_multiply(temp_h_mult_ht_matrix, prev_h_matrix, &temp_h_mult_ht_mult_h_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    for (row_index = 0; row_index < prev_h_matrix->row_count; row_index++) {
        for (column_index = 0; column_index < prev_h_matrix->column_count; column_index++) {
            /* Get (W * H(t))_ij. */
            status_code = MATRIX_get_item(temp_w_mult_h_matrix, row_index, column_index, &w_mult_h_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Get ((H(t) * H(t)^T) * H(t))_ij. */
            status_code = MATRIX_get_item(temp_h_mult_ht_mult_h_matrix, row_index, column_index, &h_mult_ht_mult_h_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Get (H(t))_ij. */
            status_code = MATRIX_get_item(prev_h_matrix, row_index, column_index, &prev_h_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Compute next (H(t+1))_ij. */
            next_h_item = prev_h_item * (1 - BETA_VALUE + BETA_VALUE * (w_mult_h_item / h_mult_ht_mult_h_item));
            status_code = MATRIX_set_item(temp_next_h_matrix, row_index, column_index, next_h_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_next_h_matrix) {
        *output_next_h_matrix = temp_next_h_matrix;
        temp_next_h_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_w_mult_h_matrix) {
        MATRIX_free(temp_w_mult_h_matrix);
        temp_w_mult_h_matrix = NULL;
    }

    if (NULL != temp_h_transpose_matrix) {
        MATRIX_free(temp_h_transpose_matrix);
        temp_h_transpose_matrix = NULL;
    }

    if (NULL != temp_h_mult_ht_matrix) {
        MATRIX_free(temp_h_mult_ht_matrix);
        temp_h_mult_ht_matrix = NULL;
    }

    if (NULL != temp_h_mult_ht_mult_h_matrix) {
        MATRIX_free(temp_h_mult_ht_mult_h_matrix);
        temp_h_mult_ht_mult_h_matrix = NULL;
    }

    if (NULL != temp_next_h_matrix) {
        MATRIX_free(temp_next_h_matrix);
        temp_next_h_matrix = NULL;
    }
    
    return status_code;
}


int SYMNMF_symnmf_algorithm(
    struct matrix *initial_h_matrix, 
    struct matrix *norm_matrix, 
    struct matrix **output_result_h_matrix
) {
    int status_code = STANDARD_ERROR_CODE;
    struct matrix *temp_prev_h_matrix = initial_h_matrix;
    struct matrix *temp_next_h_matrix = NULL;
    size_t current_iterations = 0;
    double convergence_distance = 0;

    do {
        temp_next_h_matrix = NULL;

        if (current_iterations == MAX_ITERATIONS) {
            /* Reached maximum iterations. */
            break;
        }

        /* Compute the next H matrix. */
        status_code = symnmf_update_h_matrix(temp_prev_h_matrix, norm_matrix, &temp_next_h_matrix);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        /* Compute convergence distance between prev and next H matrices. */
        convergence_distance = MATRIX_squared_frobenius_distance(temp_prev_h_matrix, temp_next_h_matrix);
        if (0 > convergence_distance) {
            goto l_cleanup;
        }

        /* Release the previous H matrix 
         * (unless it is the initial matrix, which is not owned by this function). 
         */
        if (initial_h_matrix != temp_prev_h_matrix) {
            MATRIX_free(temp_prev_h_matrix);
            temp_prev_h_matrix = NULL;   
        }

        temp_prev_h_matrix = temp_next_h_matrix;
        current_iterations++;
    } while (convergence_distance >= CONVERGENCE_THRESHOLD);

    temp_prev_h_matrix = NULL;

    /* Transfer ownership. */
    if (NULL != output_result_h_matrix) {
        *output_result_h_matrix = temp_next_h_matrix;
        temp_next_h_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    /* Ensure we are not accidentally freeing either the initial or final H matrices. */    
    if ((NULL != temp_prev_h_matrix) && (initial_h_matrix != temp_prev_h_matrix) && (temp_next_h_matrix != temp_prev_h_matrix)) {
        MATRIX_free(temp_prev_h_matrix);
        temp_prev_h_matrix = NULL;  
    }

    if (NULL != temp_next_h_matrix) {
        MATRIX_free(temp_next_h_matrix);
        temp_next_h_matrix = NULL;
    }

    return status_code;
}


int SYMNMF_calc_similarity_matrix(struct vector *datapoints_vector, size_t data_dimension, struct matrix **output_similarity_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    double item_value = 0;
    struct datapoint_item temp_first_datapoint = {0};
    struct datapoint_item temp_second_datapoint = {0};
    double datapoint_distance = 0;
    struct matrix *temp_similarity_matrix = NULL;
    size_t row_index = 0;
    size_t column_index = 0;

    /* Allocate an NxN square matrix for the similarity matrix. */
    temp_similarity_matrix = MATRIX_allocate(datapoints_vector->item_count, datapoints_vector->item_count);
    if (NULL == temp_similarity_matrix) {
        goto l_cleanup;
    }

    for (row_index = 0; row_index < temp_similarity_matrix->row_count; row_index++) {
        for (column_index = row_index; column_index < temp_similarity_matrix->column_count; column_index++) {
            item_value = 0;
            
            /* Zeros on the diagonal, otherwise we have:
             * A_ij = exp(-(||x_i - x_j||^2)/2)
             */
            if (row_index != column_index) {
                status_code = VECTOR_get_item(datapoints_vector, row_index, &temp_first_datapoint);
                if (STANDARD_SUCCESS_CODE != status_code) {
                    goto l_cleanup;
                }
                
                status_code = VECTOR_get_item(datapoints_vector, column_index, &temp_second_datapoint);
                if (STANDARD_SUCCESS_CODE != status_code) {
                    goto l_cleanup;
                }

                /* Set datapoint_distance = ||x_i - x_j||^2 */
                datapoint_distance = symnmf_squared_euclidian_distance(temp_first_datapoint.coords, temp_second_datapoint.coords, data_dimension);

                item_value = exp(-datapoint_distance/2);
            }

            /* The matrix is symmetric, so we set the transpose coordinate as well. */
            status_code = MATRIX_set_item(temp_similarity_matrix, row_index, column_index, item_value);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            status_code = MATRIX_set_item(temp_similarity_matrix, column_index, row_index, item_value);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_similarity_matrix) {
        *output_similarity_matrix = temp_similarity_matrix;
        temp_similarity_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_similarity_matrix) {
        MATRIX_free(temp_similarity_matrix);
        temp_similarity_matrix = NULL;
    }

    return status_code;
}


int SYMNMF_calc_diagonal_degree_matrix(
    struct vector *datapoints_vector,
    size_t data_dimension,
    struct matrix **output_ddg_matrix,
    struct matrix **output_similarity_matrix
) {
    int status_code = STANDARD_ERROR_CODE;
    struct matrix *temp_similarity_matrix = NULL;
    struct matrix *temp_ddg_matrix = NULL;
    double item_value = 0;
    double similarity_value = 0;
    size_t row_index = 0;
    size_t column_index = 0;
    size_t similarity_column_index = 0;

    /* Start by computing the similarity matrix. */
    status_code = SYMNMF_calc_similarity_matrix(datapoints_vector, data_dimension, &temp_similarity_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Allocate an NxN square matrix for the ddg matrix. */
    temp_ddg_matrix = MATRIX_allocate(datapoints_vector->item_count, datapoints_vector->item_count);
    if (NULL == temp_ddg_matrix) {
        goto l_cleanup;
    }

    for (row_index = 0; row_index < temp_ddg_matrix->row_count; row_index++) {
        for (column_index = row_index; column_index < temp_ddg_matrix->column_count; column_index++) {
            item_value = 0;
            
            /* Zeros off the diagonal, otherwise we have:
             * A_ii = sum(A_ij for j in range(N))
             */
            if (row_index == column_index) {
                for (similarity_column_index = 0; similarity_column_index < datapoints_vector->item_count; similarity_column_index++) {
                    status_code = MATRIX_get_item(temp_similarity_matrix, row_index, similarity_column_index, &similarity_value);
                    if (STANDARD_SUCCESS_CODE != status_code) {
                        goto l_cleanup;
                    }

                    item_value += similarity_value;
                }
            }

            /* The matrix is symmetric, so we set the transpose coordinate as well. */
            status_code = MATRIX_set_item(temp_ddg_matrix, row_index, column_index, item_value);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            status_code = MATRIX_set_item(temp_ddg_matrix, column_index, row_index, item_value);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_ddg_matrix) {
        *output_ddg_matrix = temp_ddg_matrix;
        temp_ddg_matrix = NULL;
    }

    if (NULL != output_similarity_matrix) {
        *output_similarity_matrix = temp_similarity_matrix;
        temp_similarity_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_similarity_matrix) {
        MATRIX_free(temp_similarity_matrix);
        temp_similarity_matrix = NULL;
    }

    if (NULL != temp_ddg_matrix) {
        MATRIX_free(temp_ddg_matrix);
        temp_ddg_matrix = NULL;
    }

    return status_code;
}


int SYMNMF_calc_normalized_similarity_matrix(struct vector *datapoints_vector, size_t data_dimension, struct matrix **output_norm_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    struct matrix *temp_ddg_matrix = NULL;
    struct matrix *temp_similarity_matrix = NULL;
    struct matrix *temp_multiplication_matrix = NULL;
    struct matrix *temp_norm_matrix = NULL;

    /* Start by computing the similarity and ddg matrices. */
    status_code = SYMNMF_calc_diagonal_degree_matrix(datapoints_vector, data_dimension, &temp_ddg_matrix, &temp_similarity_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Raise the ddg matrix to the (-1/2) power. */
    status_code = MATRIX_diagonal_power(temp_ddg_matrix, -0.5);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Perform first multiplication:
     * temp = (D^(-0.5)) * A 
     */
    status_code = MATRIX_multiply(temp_ddg_matrix, temp_similarity_matrix, &temp_multiplication_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Perform second multiplication:
     * W = temp * (D^(-0.5))
     */
    status_code = MATRIX_multiply(temp_multiplication_matrix, temp_ddg_matrix, &temp_norm_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Transfer ownership. */
    if (NULL != output_norm_matrix) {
        *output_norm_matrix = temp_norm_matrix;
        temp_norm_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_ddg_matrix) {
        MATRIX_free(temp_ddg_matrix);
        temp_ddg_matrix = NULL;
    }

    if (NULL != temp_similarity_matrix) {
        MATRIX_free(temp_similarity_matrix);
        temp_similarity_matrix = NULL;
    }

    if (NULL != temp_multiplication_matrix) {
        MATRIX_free(temp_multiplication_matrix);
        temp_multiplication_matrix = NULL;
    }

    if (NULL != temp_norm_matrix) {
        MATRIX_free(temp_norm_matrix);
        temp_norm_matrix = NULL;
    }

    return status_code;
}


int SYMNMF_free_datapoint_vector(struct vector *datapoints) {
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    int datapoint_index = 0;

    for (datapoint_index = datapoints->item_count - 1; datapoint_index >= 0; datapoint_index--) {
        /* Release buffer from the datastructure. */
        status_code = VECTOR_get_item(datapoints, datapoint_index, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        free(temp_datapoint_item.coords);
        temp_datapoint_item.coords = NULL;

        /* Write back reset points to the datastructure. */
        status_code = VECTOR_set_item(datapoints, datapoint_index, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        datapoints->item_count--;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    VECTOR_free(datapoints);
    return status_code;
}
