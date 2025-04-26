#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

/********************** Constants **********************/
#define DEFAULT_NUM_ITERATIONS (400)
#define CONVERGENCE_THRESHOLD (0.001)
#define MINIMUM_VALID_K_VALUE (1)
#define MINIMUM_VALID_NUM_ITERATIONS (1)
#define MAXIMUM_VALID_NUM_ITERATIONS (1000)
#define INVALID_K_VALUE_MESSAGE "Incorrect number of clusters!\n"
#define INVALID_NUM_ITERATIONS_MESSAGE "Incorrect maximum iteration!\n"
#define GENERAL_ERROR_MESSAGE "An Error Has Occurred\n"
#define DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE (32)

#define NUM_REQUIRED_ARGS (2)
#define NUM_OPTIONAL_ARGS (1)
#define K_VALUE_ARG_INDEX (1)
#define NUM_ITERATIONS_ARG_INDEX (2)

#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)
#define ERROR_EXIT_CODE (1)
#define SUCCESS_EXIT_CODE (0)

/********************** Structures **********************/
/* A dynamically increasing vector for data items of same size. 
 * Note that only adding/setting/clearing items is supported, not popping.
*/
struct vector {
    /* Allocated buffer. */
    void *data;

    /* Size of items in the vector. */
    size_t item_size;

    /* Allocated size in items. */
    size_t item_capacity;
    
    /* Occupied size in items. */
    size_t item_count;
};


/* Single datapoint item to be stored in the datapoints vector. */
struct datapoint_item {
    double *coords;
};


/* Single cluster item to be stored in the clusters vector. */
struct cluster_item {
    double *centroid;
    struct vector *datapoints;
};

/********************** Vector Operations **********************/
/**
 * Allocate a dynamically resizing vector of same-size items.
 * 
 * @param item_size The size of the items in vector, in bytes.
 * @param reserve_items The initial size to reserve for the vector.
 * @return The allocated vector or NULL if an error occurred. 
 *  The vector must be passed to free_vector() to release memory. 
 */
struct vector *allocate_vector(size_t item_size, size_t reserve_items) {
    struct vector *temp_vector = NULL;

    if ((0 == item_size) || (0 == reserve_items)) {
        return NULL;
    }

    temp_vector = malloc(sizeof(struct vector));
    if (NULL == temp_vector) {
        return NULL;
    }

    temp_vector->data = malloc(item_size * reserve_items);
    if (NULL == temp_vector->data) {
        free(temp_vector);
        return NULL;
    }

    temp_vector->item_size = item_size;
    temp_vector->item_capacity = reserve_items;
    temp_vector->item_count = 0;

    return temp_vector;
}


/**
 * Retrieve from the vector the item data at the given index. 
 * 
 * @param vector The vector object.
 * @param item_index The index of the item.
 * @param item_buffer Buffer to which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int get_vector_item(struct vector *vector, size_t item_index, void *item_buffer) {
    int status_code = STANDARD_ERROR_CODE;
    
    if (vector->item_count <= item_index) {
        /* Index out of range of initialized items. */
        goto l_cleanup;
    } else {
        memcpy(item_buffer, ((unsigned char *) vector->data) + (item_index * vector->item_size), vector->item_size);
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Set the item data in the vector at the given index. 
 * 
 * @param vector The vector object.
 * @param item_index The index of the item.
 * @param item_buffer Buffer from which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int set_vector_item(struct vector *vector, size_t item_index, void *item_data) {
    int status_code = STANDARD_ERROR_CODE;
    
    if (vector->item_count <= item_index) {
        /* Index out of range of initialized items. */
        goto l_cleanup;
    } else {
        memcpy(((unsigned char *) vector->data) + (item_index * vector->item_size), item_data, vector->item_size);
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Append an item with given data to the vector.
 * If no capacity is left in the vector, more space will be allocated. 
 * 
 * @param vector The vector object.
 * @param item_buffer Buffer from which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int add_vector_item(struct vector *vector, void *item_data) {
    int status_code = STANDARD_ERROR_CODE;

    if (vector->item_count == vector->item_capacity) {
        vector->data = realloc(vector->data, vector->item_size * (vector->item_capacity * 2));
        if (NULL == vector->data) {
            goto l_cleanup;
        }

        vector->item_capacity *= 2;
    } 

    memcpy(((unsigned char *) vector->data) + (vector->item_count * vector->item_size), item_data, vector->item_size);
    vector->item_count++;

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Release the vector object and its memory in use.
 * Note that if the items theselves point to allocated memory, their allocations will not be released. 
 * 
 * @param vector The vector object.
 */
void free_vector(struct vector *vector) {
    if (NULL == vector) {
        return;
    }

    free(vector->data);
    free(vector);
}


/********************** K-Means Algorithm **********************/
/**
 * Parse the first line of the input file to conclude the dimension of the coordinates.
 * We assume the line is valid, with comma-delimited floats, and an equal amount of coordinates for all datapoints. 
 * 
 * @param datapoint_buffer The first line of datapoint input.
 * @param datapoint_length The length of the first line.
 * @return The amount of commas seen in the first line.
 */
size_t parse_data_dimension(const char *datapoint_buffer, size_t datapoint_length) {
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
 * Parse the given input line and append a new datapoint into the datapoints vector.
 * We assume all input lines are of the same length as the first line, i.e. same dimension.
 * 
 * @param datapoint_buffer The datapoint input line to process. Must be at least the length of the first input line.
 * @param data_dimension The dimension of the datapoints (i.e. number of coordinates per point)
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int add_datapoint_to_vector(
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
    status_code = add_vector_item(datapoints_vector, &temp_datapoint_item);
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
 * Release the clusters vector and all its associated allocations.
 * 
 * @param clusters The clusters vector to free.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int free_cluster_vector(struct vector *clusters) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    int i = 0;

    for (i = clusters->item_count - 1; i >= 0; i--) {
        /* Release buffer from the datastructure. */
        status_code = get_vector_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        free(temp_cluster_item.centroid);
        temp_cluster_item.centroid = NULL;

        /* Release cluster vector but not the buffers used for their datapoints,
         * because points are owned by the datapoints vector.
         */
        free_vector(temp_cluster_item.datapoints);
        temp_cluster_item.datapoints = NULL;

        /* Write back reset points to the datastructure. */
        status_code = set_vector_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        clusters->item_count--;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free_vector(clusters);
    return status_code;
}


/**
 * Release the datapoints vector and all its associated allocations.
 * 
 * @param datapoints The datapoints vector to free.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int free_datapoint_vector(struct vector *datapoints) {
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    int i = 0;

    for (i = datapoints->item_count - 1; i >= 0; i--) {
        /* Release buffer from the datastructure. */
        status_code = get_vector_item(datapoints, i, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        free(temp_datapoint_item.coords);
        temp_datapoint_item.coords = NULL;

        /* Write back reset points to the datastructure. */
        status_code = set_vector_item(datapoints, i, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        datapoints->item_count--;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free_vector(datapoints);
    return status_code;
}


/**
 * Set the initial values of the clusters: Empty cluster vectors, and centroids set to the first K datapoints.
 * 
 * @param clusters The clusters vector to initialize.
 * @param datapoints The datapoints vector loaded with input datapoints.
 * @param k_value The K parameters used for the K-means algorithm.
 * @param dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int initialize_clusters(struct vector *clusters, struct vector *datapoints, size_t k_value, size_t dimension) {
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    struct cluster_item temp_cluster_item = {0};
    double *temp_coords = NULL;
    struct vector *temp_cluster_vector = NULL;
    size_t i = 0;

    /* Set first K datapoints as clusters. */
    for (i = 0; i < k_value; i++) {
        /* Allocate buffer for centroid coords. */
        temp_coords = malloc(sizeof(double) * dimension);
        if (NULL == temp_coords) {
            goto l_cleanup;
        }
        temp_cluster_item.centroid = temp_coords;

        /* Read coord data from datapoints. */
        status_code = get_vector_item(datapoints, i, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }
        memcpy(temp_cluster_item.centroid, temp_datapoint_item.coords, sizeof(double) * dimension);

        /* Allocate vector for cluster. */
        temp_cluster_vector = allocate_vector(sizeof(struct datapoint_item), DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE);
        if (NULL == temp_cluster_vector) {
            goto l_cleanup;
        }
        temp_cluster_item.datapoints = temp_cluster_vector;

        /* Insert the centroid into the vector. */
        status_code = add_vector_item(clusters, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        temp_coords = NULL;
        temp_cluster_vector = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(temp_coords);
    free_vector(temp_cluster_vector);
    return status_code;
}


/**
 * Calculate the Euclidian distance of two given datapoints (of the same dimension).
 * 
 * @param first_point The first datapoint in the calculation.
 * @param second_point The second datapoint in the calculation.
 * @param dimension The dimension of the datapoints.
 * @return The computed distance.
 */
double calculate_euclidian_distance(double *first_point, double *second_point, size_t dimension) {
    double temp_distance = 0;
    size_t i = 0;

    for (i = 0; i < dimension; i++) {
        temp_distance += pow(first_point[i] - second_point[i], 2);
    }

    return sqrt(temp_distance);
}


/**
 * Iterate over the clusters to find the one with the closest distance to a given point.
 * 
 * @param clusters The clusters vector.
 * @param datapoint The datapoint to check against.
 * @param data_dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int find_closest_cluster(
    struct vector *clusters, 
    struct datapoint_item *datapoint, 
    size_t data_dimension
) {
    int status_code = STANDARD_ERROR_CODE;
    double minimum_distance = 0;
    double temp_distance = 0;
    int minimum_cluster = -1;
    struct cluster_item temp_cluster_item = {0};
    size_t i = 0;

    for (i = 0; i < clusters->item_count; i++) {
        status_code = get_vector_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        temp_distance = calculate_euclidian_distance(datapoint->coords, temp_cluster_item.centroid, data_dimension);
        if ((temp_distance < minimum_distance) || (0 == i)) {
            minimum_distance = temp_distance;
            minimum_cluster = i;
        }
    }

    if (minimum_cluster < 0) {
        goto l_cleanup;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (STANDARD_SUCCESS_CODE == status_code) {
        return minimum_cluster;
    } else {
        return STANDARD_ERROR_CODE;
    }
}


/**
 * Insert a given datapoint to a cluster's datapoint list.
 * 
 * @param clusters The clusters vector.
 * @param datapoint The datapoint to insert.
 * @param cluster_index The cluster index to insert into.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int assign_datapoint_to_cluster(
    struct vector *clusters, 
    struct datapoint_item *datapoint,
    size_t cluster_index
) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};

    status_code = get_vector_item(clusters, cluster_index, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    status_code = add_vector_item(temp_cluster_item.datapoints, datapoint);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Clear all datapoints from all clusters' datapoint lists.
 * 
 * @param clusters The clusters vector.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int remove_cluster_assignments(struct vector *clusters) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    size_t i = 0;

    for (i = 0; i < clusters->item_count; i++) {
        status_code = get_vector_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        /* Reset the count without removing underlying memory, datapoints are owned by datapoint vector. */
        temp_cluster_item.datapoints->item_count = 0;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Recalculate the centroid of the datapoints assigned to a given cluster.
 * 
 * @param clusters The clusters vector.
 * @param cluster_index The cluster index to update.
 * @param dimension The dimension of the datapoints.
 * @param output_diff An output parameter to hold the delta in distance between the previous and updated centroids. 
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int update_cluster_centroid(struct vector *clusters, size_t cluster_index, size_t dimension, double *output_diff) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    struct datapoint_item temp_datapoint_item = {0};
    size_t dimension_index = 0;
    size_t item_index = 0;
    double *previous_centroid = NULL;

    status_code = get_vector_item(clusters, cluster_index, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Save the previous centroid so we can check the diff. */
    previous_centroid = malloc(sizeof(double) * dimension);
    if (NULL == previous_centroid) {
        goto l_cleanup;
    }
    memcpy(previous_centroid, temp_cluster_item.centroid, sizeof(double) * dimension);

    /* Calculate centroid of cluster using datapoints in the cluster. */
    for (dimension_index = 0; dimension_index < dimension; dimension_index++) {
        temp_cluster_item.centroid[dimension_index] = 0;

        for (item_index = 0; item_index < temp_cluster_item.datapoints->item_count; item_index++) {
            status_code = get_vector_item(temp_cluster_item.datapoints, item_index, &temp_datapoint_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            temp_cluster_item.centroid[dimension_index] += temp_datapoint_item.coords[dimension_index];
        }

        temp_cluster_item.centroid[dimension_index] /= temp_cluster_item.datapoints->item_count;
    }
    
    /* Write back results to the cluster item. */
    status_code = set_vector_item(clusters, cluster_index, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    
    /* Calculate and report back the diff. */
    *output_diff = calculate_euclidian_distance(previous_centroid, temp_cluster_item.centroid, dimension);

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(previous_centroid);
    return status_code;
}


/**
 * Execute the K-means algorithm with the given configurations.
 * 
 * @param clusters The empty clusters vector.
 * @param datapoints The preloaded datapoints vector.
 * @param k_value The K parameter for the algorithm (cluster count).
 * @param num_iterations The iter parameter for the algorithm (maximum iterations).
 * @param data_dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int perform_kmeans_clustering(
    struct vector *clusters,
    struct vector *datapoints,
    size_t k_value,
    size_t num_iterations,
    size_t data_dimension
) {
    int status_code = STANDARD_ERROR_CODE;
    size_t current_iterations = 0;
    double maximum_convergence_distance = 0;
    int closest_cluster = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    size_t datapoint_index = 0;
    size_t cluster_index = 0;
    double update_centroid_diff = 0;

    /* Initialize clusters as first k datapoints. */
    status_code = initialize_clusters(clusters, datapoints, k_value, data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    do {
        maximum_convergence_distance = 0;

        if (current_iterations == num_iterations) {
            /* Reached maximum iterations. */
            break;
        }
        
        /* Remove all datapoint assignments. */
        status_code = remove_cluster_assignments(clusters);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        /* Assign every datapoint to the closest cluster. */
        for (datapoint_index = 0; datapoint_index < datapoints->item_count; datapoint_index++) {
            status_code = get_vector_item(datapoints, datapoint_index, &temp_datapoint_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            closest_cluster = find_closest_cluster(clusters, &temp_datapoint_item, data_dimension);
            if (STANDARD_ERROR_CODE == closest_cluster) {
                goto l_cleanup;
            }

            status_code = assign_datapoint_to_cluster(clusters, &temp_datapoint_item, closest_cluster);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }

        /* Update centroids. */
        for (cluster_index = 0; cluster_index < clusters->item_count; cluster_index++) {
            status_code = update_cluster_centroid(clusters, cluster_index, data_dimension, &update_centroid_diff);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Maintain the largest centroid diff to check convergence. */
            if (update_centroid_diff > maximum_convergence_distance) {
                maximum_convergence_distance = update_centroid_diff;
            }
        }

        current_iterations++;
    } while (maximum_convergence_distance >= CONVERGENCE_THRESHOLD);

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


/**
 * Print the output of the K-means algorithm using the implementation's datastructures.
 * 
 * @param clusters The clusters vector with the final clusters.
 * @param dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int print_kmeans_results(struct vector *clusters, size_t dimension) {
    int status_code = STANDARD_ERROR_CODE;
    size_t cluster_index = 0;
    size_t dimension_index = 0;
    struct cluster_item temp_cluster_item = {0};

    for (cluster_index = 0; cluster_index < clusters->item_count; cluster_index++) {
        status_code = get_vector_item(clusters, cluster_index, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        for (dimension_index = 0; dimension_index < dimension; dimension_index++) {
            /* Print coordinates of centroid. */
            printf("%.4f", temp_cluster_item.centroid[dimension_index]);

            /* Print delimiter if not the last coord. */
            if (dimension_index < dimension - 1) {
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
 * Program entry point. Load input into datastructures, execute K-means algorithm and print out results.
 * Usage: ./kmeans <k-param> [iter-param]
 * Input datapoints received via stdin. Input is assumed to be valid.
 * 
 * @return The status of the program. SUCCESS_EXIT_CODE for success, ERROR_EXIT_CODE for error.
 */
int main(int argc, char *argv[]) {
    size_t k_value = 0;
    size_t num_iterations = DEFAULT_NUM_ITERATIONS;
    size_t data_dimension = 0;
    char *line_buffer = NULL;
    size_t line_length = 0;
    ssize_t num_bytes_read = STANDARD_ERROR_CODE;
    struct vector *datapoints_vector = NULL;
    struct vector *clusters_vector = NULL; 
    int exit_code = ERROR_EXIT_CODE;
    int status_code = STANDARD_ERROR_CODE;

    /* Parse required command line arguments. */
    if (NUM_REQUIRED_ARGS > argc) {
        printf(GENERAL_ERROR_MESSAGE);
        goto l_cleanup;
    }

    k_value = strtoul(argv[K_VALUE_ARG_INDEX], NULL, 10);
    if (MINIMUM_VALID_K_VALUE > k_value) {
        printf(INVALID_K_VALUE_MESSAGE);
        goto l_cleanup;
    }

    /* Check for optional command line arguments. */
    if (NUM_REQUIRED_ARGS < argc) {
        if ((NUM_REQUIRED_ARGS + NUM_OPTIONAL_ARGS) < argc) {
            printf(GENERAL_ERROR_MESSAGE);
            goto l_cleanup;
        }

        num_iterations = strtoul(argv[NUM_ITERATIONS_ARG_INDEX], NULL, 10);
        if ((MINIMUM_VALID_NUM_ITERATIONS > num_iterations) || (MAXIMUM_VALID_NUM_ITERATIONS < num_iterations)) {
            printf(INVALID_NUM_ITERATIONS_MESSAGE);
            goto l_cleanup;
        } 
    } 

    /* Read first line of input to determine datapoint dimensions. */
    num_bytes_read = getline(&line_buffer, &line_length, stdin);
    if (0 >= num_bytes_read) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }
    data_dimension = parse_data_dimension(line_buffer, num_bytes_read);

    /* Initialize datastructures from parsed dimension. */
    datapoints_vector = allocate_vector(sizeof(struct datapoint_item), DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE);
    if (NULL == datapoints_vector) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }    

    clusters_vector = allocate_vector(sizeof(struct cluster_item), k_value);
    if (NULL == clusters_vector) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }  

    /* Parse and load all datapoints into vector. */
    do {
        status_code = add_datapoint_to_vector(line_buffer, data_dimension, datapoints_vector);
        if (STANDARD_SUCCESS_CODE != status_code) {
            printf(GENERAL_ERROR_MESSAGE);  
            goto l_cleanup;
        }

        /* Reset and check errno after call to see if error occurred or EOF reached.
         * If there is no error, we have reached the end of the input.
         */
        errno = 0;
        num_bytes_read = getline(&line_buffer, &line_length, stdin);
        if ((STANDARD_ERROR_CODE == num_bytes_read) && (0 != errno)) {
            printf(GENERAL_ERROR_MESSAGE);  
            goto l_cleanup;
        }
    } while (num_bytes_read > 0);

    /* Count datapoints and validate k-value. */
    if (datapoints_vector->item_count <= k_value) {
        printf(INVALID_K_VALUE_MESSAGE);
        goto l_cleanup;
    }

    /* Execute k-means algorithm. */
    status_code = perform_kmeans_clustering(
        clusters_vector,
        datapoints_vector,
        k_value,
        num_iterations, 
        data_dimension
    );
    if (STANDARD_SUCCESS_CODE != status_code) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    }

    /* Print results. */
    status_code = print_kmeans_results(clusters_vector, data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        printf(GENERAL_ERROR_MESSAGE);  
        goto l_cleanup;
    } 

    exit_code = SUCCESS_EXIT_CODE;
l_cleanup:
    if (NULL != clusters_vector) {
        free_cluster_vector(clusters_vector);
    }
    if (NULL != datapoints_vector) {
        free_datapoint_vector(datapoints_vector);
    }
    free(line_buffer);
    return exit_code;
}
