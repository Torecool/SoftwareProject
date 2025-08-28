#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "symnmf.h"
#include "vector.h"
#include "matrix.h"

/********************** Constants *****************************/
#define GENERAL_ERROR_MESSAGE "An Error Has Occurred\n"
#define DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE (32)

/* Num required args including program name. */
#define NUM_REQUIRED_ARGS (3)
#define GOAL_ARG_INDEX (1)
#define FILENAME_ARG_INDEX (2)

#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)
#define ERROR_EXIT_CODE (1)
#define SUCCESS_EXIT_CODE (0)


/********************** Program Entrypoint ********************/
/**
 * Parse the first line of the input file to conclude the dimension of the coordinates.
 * We assume the line is valid, with comma-delimited floats, and an equal amount of coordinates for all datapoints. 
 * 
 * @param datapoint_buffer The first line of datapoint input.
 * @param datapoint_length The length of the first line.
 * @return The amount of commas seen in the first line.
 */
static size_t symnmf_parse_data_dimension(const char *datapoint_buffer, size_t datapoint_length) {
    /* Minimum is 1, one comma delimiter separates two coords. */
    size_t count = 1;
    size_t i = 0;

    for (i = 0; i < datapoint_length; i++) {
        if (datapoint_buffer[i] == ',') {
            count++;
        }
    }

    return count;
}


/**
 * Print the output of the SymNMF algorithm using the implementation's datastructures.
 * 
 * @param result_matrix The matrix object with the output data.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int symnmf_print_results(struct matrix *result_matrix) {
    int status_code = STANDARD_ERROR_CODE;
    double item_data = 0;
    size_t row_index = 0;
    size_t column_index = 0;

    for (row_index = 0; row_index < result_matrix->row_count; row_index++) {
        for (column_index = 0; column_index < result_matrix->column_count; column_index++) {
            status_code = MATRIX_get_item(result_matrix, row_index, column_index, &item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Print matrix item. */
            printf("%.4f", item_data);

            /* Print delimiter if not the last column. */
            if (column_index < result_matrix->column_count - 1) {
                printf(",");
            }
        }

        printf("\n");
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Parse the given input line and append a new datapoint into the datapoints vector.
 * We assume all input lines are of the same length as the first line, i.e. same dimension.
 * 
 * @param datapoint_buffer The datapoint input line to process. Must be at least the length of the first input line.
 * @param data_dimension The dimension of the datapoints (i.e. number of coordinates per point)
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int symnmf_add_datapoint_to_vector(
    const char *datapoint_buffer,
    size_t data_dimension,
    struct vector *datapoints_vector
) {
    int status_code = STANDARD_ERROR_CODE;
    const char *buffer_cursor = datapoint_buffer;
    char *temp_cursor = NULL;
    struct datapoint_item temp_datapoint_item = {0};
    size_t i = 0;

    double *temp_data = malloc(sizeof(double) * data_dimension);
    if (NULL == temp_data) {
        goto l_cleanup;
    }

    for (i = 0; i < data_dimension; i++) {
        temp_data[i] = strtod(buffer_cursor, &temp_cursor);

        /* Skip the comma delimiter. */
        buffer_cursor = temp_cursor + 1;
    }

    /* Insert the datapoint into the datastructure. */
    temp_datapoint_item.coords = temp_data;
    status_code = VECTOR_add_item(datapoints_vector, &temp_datapoint_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    temp_data = NULL;

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(temp_data);
    return status_code;
}


/**
 * Parse the datapoints input file and load into a data structure. 
 * 
 * @param filename The file path of the datapoints file.
 * @param output_datapoints_vector An optional output parameter to contain the data structure.
 *  Must later be freed using SYMNMF_free_datapoint_vector.
 * @param output_data_dimension An optional output parameter to contain the dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int symnmf_parse_datapoints_file(
    const char *filename, 
    struct vector **output_datapoints_vector, 
    size_t *output_data_dimension
) {
    int status_code = STANDARD_ERROR_CODE;
    FILE *input_file = NULL;
    char *line_buffer = NULL;
    size_t line_length = 0;
    ssize_t num_bytes_read = STANDARD_ERROR_CODE;
    size_t data_dimension = 0;
    struct vector *temp_datapoints_vector = NULL;

    /* Allocate a data structure for the datapoints. */
    temp_datapoints_vector = VECTOR_allocate(sizeof(struct datapoint_item), DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE);
    if (NULL == temp_datapoints_vector) {
        goto l_cleanup;
    }

    /* Try to open input file. */
    input_file = fopen(filename, "r");
    if (NULL == input_file) {
        goto l_cleanup;
    }

    /* Read first line of input to determine datapoint dimensions. */
    num_bytes_read = getline(&line_buffer, &line_length, input_file);
    if (0 >= num_bytes_read) {
        goto l_cleanup;
    }
    data_dimension = symnmf_parse_data_dimension(line_buffer, num_bytes_read);   

    /* Parse and load all datapoints into vector. */
    do {
        status_code = symnmf_add_datapoint_to_vector(line_buffer, data_dimension, temp_datapoints_vector);
        if (STANDARD_SUCCESS_CODE != status_code) {
            printf(GENERAL_ERROR_MESSAGE);  
            goto l_cleanup;
        }

        /* Reset and check errno after call to see if error occurred or EOF reached.
         * If there is no error, we have reached the end of the input.
         */
        errno = 0;
        num_bytes_read = getline(&line_buffer, &line_length, input_file);
        if ((STANDARD_ERROR_CODE == num_bytes_read) && (0 != errno)) {
            printf(GENERAL_ERROR_MESSAGE);  
            goto l_cleanup;
        }
    } while (num_bytes_read > 0);

    /* Transfer ownership. */
    if (NULL != output_datapoints_vector) {
        *output_datapoints_vector = temp_datapoints_vector;
        temp_datapoints_vector = NULL;
    }

    if (NULL != output_data_dimension) {
        *output_data_dimension = data_dimension;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != input_file) {
        fclose(input_file);
    }
    
    free(line_buffer);

    if (NULL != temp_datapoints_vector) {
        SYMNMF_free_datapoint_vector(temp_datapoints_vector);
        temp_datapoints_vector = NULL;
    }

    return status_code;
}


/**
 * Execute the goal matching the given name.
 * 
 * @param goal_name The name of the goal to execute.
 * @param datapoints_vector The datapoints structure.
 * @param data_dimension The dimension of the datapoints.
 * @param output_matrix An optional output parameter to contain matrix which is the result of the operation.
 *  Must later be freed via MATRIX_free().
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int execute_program_goal(
    const char *goal_name, 
    struct vector *datapoints_vector, 
    size_t data_dimension,
    struct matrix **output_matrix
) {
    int status_code = STANDARD_ERROR_CODE;

    if (0 == strcmp("sym", goal_name)) {
        status_code = SYMNMF_calc_similarity_matrix(datapoints_vector, data_dimension, output_matrix);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }
    } else if (0 == strcmp("ddg", goal_name)) {
        status_code = SYMNMF_calc_diagonal_degree_matrix(datapoints_vector, data_dimension, output_matrix, NULL);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }
    } else if (0 == strcmp("norm", goal_name)) {
        status_code = SYMNMF_calc_normalized_similarity_matrix(datapoints_vector, data_dimension, output_matrix);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }
    } else {
        /* Invalid goal. */
        goto l_cleanup;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Program entry point. Load input into datastructures, execute SymNMF algorithm and print out results.
 * Usage: ./symnmf <goal-param> <filename-param>
 * Input is assumed to be valid.
 * 
 * @return The status of the program. SUCCESS_EXIT_CODE for success, ERROR_EXIT_CODE for error.
 */
int main(int argc, char *argv[]) {
    const char *goal_param = NULL;
    const char *filename_param = NULL;
    size_t data_dimension = 0;
    struct vector *datapoints_vector = NULL;
    struct matrix *result_matrix = NULL;
    int exit_code = ERROR_EXIT_CODE;
    int status_code = STANDARD_ERROR_CODE;

    /* Parse required command line arguments. */
    if (NUM_REQUIRED_ARGS != argc) {
        printf(GENERAL_ERROR_MESSAGE);
        goto l_cleanup;
    }

    goal_param = argv[GOAL_ARG_INDEX];
    filename_param = argv[FILENAME_ARG_INDEX];

    /* Load input data. */
    status_code = symnmf_parse_datapoints_file(filename_param, &datapoints_vector, &data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }

    /* Execute program goal. */
    status_code = execute_program_goal(goal_param, datapoints_vector, data_dimension, &result_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }

    /* Print results. */
    status_code = symnmf_print_results(result_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    } 

    exit_code = SUCCESS_EXIT_CODE;
l_cleanup:
    if (NULL != datapoints_vector) {
        SYMNMF_free_datapoint_vector(datapoints_vector);
    }

    if (NULL != result_matrix) {
        MATRIX_free(result_matrix);
    }
    return exit_code;
}
