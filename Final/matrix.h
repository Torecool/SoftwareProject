#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdlib.h>

/********************** Constants ******************************/
#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)

/********************** Structures *****************************/
/* A fixed-size matrix of doubles. */
struct matrix {
    /* Allocated buffer. */
    double *data;

    /* Number of rows in the matrix. */
    size_t row_count;
    
    /* Number of columns in the matrix. */
    size_t column_count;
};

/********************** Functions ******************************/
/**
 * Allocate a fixed-size matrix of doubles.
 * 
 * @param row_count Number of rows in the matrix.
 * @param column_count Number of cloumns in the matrix.
 * @return The allocated matrix or NULL if an error occurred. 
 *  The matrix must be passed to MATRIX_free() to release memory. 
 */
struct matrix *MATRIX_allocate(size_t row_count, size_t column_count);

/**
 * Retrieve from the matrix the item data at the given row and column. 
 * 
 * @param matrix The matrix object.
 * @param row_index The row of the item.
 * @param column_index The column of the item.
 * @param output_item Output variable to which the item data will be copied.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int MATRIX_get_item(struct matrix *matrix, size_t row_index, size_t column_index, double *output_item);

/**
 * Set the item data in the matrix at the given row and column. 
 * 
 * @param matrix The matrix object.
 * @param row_index The row of the item.
 * @param column_index The column of the item.
 * @param output_item Output variable from which the item data will be copied.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int MATRIX_set_item(struct matrix *matrix, size_t row_index, size_t column_index, double item_data);

/**
 * Raise the items on the diagonal to a given power of a given exponent.
 * The matrix must be square and diagonal.
 * The multiplication is performed in-place.
 * 
 * @param matrix The matrix object.
 * @param exponent The exponent to apply to the diagonal.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int MATRIX_diagonal_power(struct matrix *matrix, double exponent);

/**
 * Multiply two given matrices and output the resulting matrix.
 * The column count of the first matrix must match the row count of the second matrix.
 * 
 * @param first_matrix The first matrix object in the multiplication.
 * @param second_matrix The second matrix object in the multiplication.
 * @param output_matrix Output variable to which the resulting matrix will be transferred.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 *  The output matrix must be passed to MATRIX_free() to release memory. 
 */
int MATRIX_multiply(struct matrix *first_matrix, struct matrix *second_matrix, struct matrix **output_matrix);

/**
 * Transpose the given matrix and output the resulting matrix.
 * 
 * @param matrix The matrix object in the calculation.
 * @param output_matrix Output variable to which the resulting matrix will be transferred.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 *  The output matrix must be passed to MATRIX_free() to release memory. 
 */
int MATRIX_transpose(struct matrix *matrix, struct matrix **output_matrix);

/**
 * Calculate the sqaured Frobenius distance between two given matrices.
 * The matrices must have the same dimensions.
 * 
 * @param first_matrix The first matrix object in the calculation.
 * @param second_matrix The second matrix object in the calculation.
 * @return The computed norm, or a negative value if an error occurred.
 */
double MATRIX_squared_frobenius_distance(struct matrix *first_matrix, struct matrix *second_matrix);

/**
 * Release the matrix object and its memory in use.
 * 
 * @param matrix The matrix object.
 */
void MATRIX_free(struct matrix *matrix);

#endif /* __MATRIX_H__ */
