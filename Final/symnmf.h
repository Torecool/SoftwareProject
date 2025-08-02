#ifndef __SYMNMF_H__
#define __SYMNMF_H__

#include "vector.h"
#include "matrix.h"

/********************** Structures ****************************/
/* Single datapoint item to be stored in the datapoints vector. */
struct datapoint_item {
    double *coords;
};

/********************** Functions *****************************/
/**
 * Execute the SymNMF algorithm to compute the final H matrix.
 * 
 * @param initial_h_matrix The initial H matrix to be used in the algorithm.
 * @param norm_matrix The normalized similarity W matrix to be used in the algorithm.
 * @param output_result_h_matrix An optional output variable to which the resulting H matrix will be transferred.
 *  Must be passed to MATRIX_free() to release memory. 
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int SYMNMF_symnmf_algorithm(
    struct matrix *initial_h_matrix, 
    struct matrix *norm_matrix, 
    struct matrix **output_result_h_matrix
);

/**
 * Calculate and output the similarity matrix.
 * 
 * @param datapoints_vector The input datapoints vector.
 * @param data_dimension The dimension of the datapoints.
 * @param output_similarity_matrix An optional output variable to which the resulting similarity matrix will be transferred.
 *  Must be passed to MATRIX_free() to release memory. 
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int SYMNMF_calc_similarity_matrix(
    struct vector *datapoints_vector, 
    size_t data_dimension,
    struct matrix **output_similarity_matrix
);

/**
 * Calculate and output the diagonal degree matrix.
 * 
 * @param datapoints_vector The input datapoints vector.
 * @param data_dimension The dimension of the datapoints.
 * @param output_ddg_matrix An optional output variable to which the resulting ddg matrix will be transferred.
 *  Must be passed to MATRIX_free() to release memory. 
 * @param output_similarity_matrix An optinal output variable to which the resulting similarity matrix will be transferred.
 *  Must be passed to MATRIX_free() to release memory. 
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int SYMNMF_calc_diagonal_degree_matrix(
    struct vector *datapoints_vector,
    size_t data_dimension,
    struct matrix **output_ddg_matrix,
    struct matrix **output_similarity_matrix
);

/**
 * Calculate and output the normalized similarity matrix.
 * 
 * @param datapoints_vector The input datapoints vector.
 * @param data_dimension The dimension of the datapoints.
 * @param output_norm_matrix An optional output variable to which the resulting norm matrix will be transferred.
 *  Must be passed to MATRIX_free() to release memory. 
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int SYMNMF_calc_normalized_similarity_matrix(
    struct vector *datapoints_vector, 
    size_t data_dimension,
    struct matrix **output_norm_matrix
);

/**
 * Release the datapoints vector and all its associated allocations.
 * 
 * @param datapoints The datapoints vector to free.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int SYMNMF_free_datapoint_vector(struct vector *datapoints);

#endif /* __KMEANS_H__ */
