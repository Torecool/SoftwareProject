#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "vector.h"
#include "kmeans.h"

/********************** Functions *****************************/
/**
 * Calculate the Euclidian distance of two given datapoints (of the same dimension).
 * 
 * @param first_point The first datapoint in the calculation.
 * @param second_point The second datapoint in the calculation.
 * @param dimension The dimension of the datapoints.
 * @return The computed distance.
 */
static double kmeans_calculate_euclidian_distance(double *first_point, double *second_point, size_t dimension) {
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
static int kmeans_find_closest_cluster(
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
        status_code = VECTOR_get_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        temp_distance = kmeans_calculate_euclidian_distance(datapoint->coords, temp_cluster_item.centroid, data_dimension);
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
 * Clear all datapoints from all clusters' datapoint lists.
 * 
 * @param clusters The clusters vector.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int kmeans_remove_cluster_assignments(struct vector *clusters) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    size_t i = 0;

    for (i = 0; i < clusters->item_count; i++) {
        status_code = VECTOR_get_item(clusters, i, &temp_cluster_item);
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
static int kmeans_update_cluster_centroid(struct vector *clusters, size_t cluster_index, size_t dimension, double *output_diff) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    struct datapoint_item temp_datapoint_item = {0};
    size_t dimension_index = 0;
    size_t item_index = 0;
    double *previous_centroid = NULL;

    status_code = VECTOR_get_item(clusters, cluster_index, &temp_cluster_item);
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
            status_code = VECTOR_get_item(temp_cluster_item.datapoints, item_index, &temp_datapoint_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            temp_cluster_item.centroid[dimension_index] += temp_datapoint_item.coords[dimension_index];
        }

        temp_cluster_item.centroid[dimension_index] /= temp_cluster_item.datapoints->item_count;
    }
    
    /* Write back results to the cluster item. */
    status_code = VECTOR_set_item(clusters, cluster_index, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    
    /* Calculate and report back the diff. */
    *output_diff = kmeans_calculate_euclidian_distance(previous_centroid, temp_cluster_item.centroid, dimension);

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(previous_centroid);
    return status_code;
}


/**
 * Insert a given datapoint to a cluster's datapoint list.
 * 
 * @param clusters The clusters vector.
 * @param datapoint The datapoint to insert.
 * @param cluster_index The cluster index to insert into.
 * @param data_dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
static int kmeans_assign_datapoint_to_cluster(
    struct vector *clusters, 
    struct datapoint_item *datapoint,
    size_t cluster_index,
    size_t data_dimension
) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    
    status_code = VECTOR_get_item(clusters, cluster_index, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    status_code = VECTOR_add_item(temp_cluster_item.datapoints, datapoint);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int KMEANS_create_cluster_item(struct vector *clusters, struct datapoint_item *initial_centroid, size_t data_dimension) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    struct vector *temp_cluster_vector = NULL;
    void *temp_coords = NULL;

    /* Allocate memory for cluster item data. */
    temp_coords = calloc(data_dimension, sizeof(double));
    if (NULL == temp_coords) {
        goto l_cleanup;
    }
    temp_cluster_item.centroid = temp_coords;

    memcpy(temp_cluster_item.centroid, initial_centroid->coords, sizeof(double) * data_dimension);

    /* Allocate vector for cluster datapoints. */
    temp_cluster_vector = VECTOR_allocate(sizeof(struct datapoint_item), DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE);
    if (NULL == temp_cluster_vector) {
        goto l_cleanup;
    }
    temp_cluster_item.datapoints = temp_cluster_vector;
    temp_coords = NULL;

    status_code = VECTOR_add_item(clusters, &temp_cluster_item);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }
    temp_cluster_vector = NULL;

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(temp_coords);

    if (NULL != temp_cluster_vector) {
        KMEANS_free_cluster_vector(temp_cluster_vector);
    }

    return status_code;
}


int KMEANS_free_cluster_vector(struct vector *clusters) {
    int status_code = STANDARD_ERROR_CODE;
    struct cluster_item temp_cluster_item = {0};
    int i = 0;

    if (NULL == clusters) {
        status_code = STANDARD_SUCCESS_CODE;
        return status_code;
    }

    for (i = clusters->item_count - 1; i >= 0; i--) {
        /* Release buffer from the datastructure. */
        status_code = VECTOR_get_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        free(temp_cluster_item.centroid);
        temp_cluster_item.centroid = NULL;

        /* Release cluster vector but not the buffers used for their datapoints,
         * because points are owned by the datapoints vector.
         */
        VECTOR_free(temp_cluster_item.datapoints);
        temp_cluster_item.datapoints = NULL;

        /* Write back reset points to the datastructure. */
        status_code = VECTOR_set_item(clusters, i, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        clusters->item_count--;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    VECTOR_free(clusters);
    return status_code;
}


int KMEANS_free_datapoint_vector(struct vector *datapoints) {
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    int i = 0;

    for (i = datapoints->item_count - 1; i >= 0; i--) {
        /* Release buffer from the datastructure. */
        status_code = VECTOR_get_item(datapoints, i, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        free(temp_datapoint_item.coords);
        temp_datapoint_item.coords = NULL;

        /* Write back reset points to the datastructure. */
        status_code = VECTOR_set_item(datapoints, i, &temp_datapoint_item);
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


int KMEANS_perform_kmeans_clustering(
    struct vector *clusters,
    struct vector *datapoints,
    size_t k_value,
    size_t num_iterations,
    size_t data_dimension,
    double convergence_threshold
) {
    int status_code = STANDARD_ERROR_CODE;
    size_t current_iteration = 0;
    double maximum_convergence_distance = 0;
    int closest_cluster = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    size_t datapoint_index = 0;
    size_t cluster_index = 0;
    double update_centroid_diff = 0;

    do {
        maximum_convergence_distance = 0;

        if (current_iteration == num_iterations) {
            /* Reached maximum iterations. */
            break;
        }
        
        /* Remove all datapoint assignments. */
        status_code = kmeans_remove_cluster_assignments(clusters);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        /* Assign every datapoint to the closest cluster. */
        for (datapoint_index = 0; datapoint_index < datapoints->item_count; datapoint_index++) {
            status_code = VECTOR_get_item(datapoints, datapoint_index, &temp_datapoint_item);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            closest_cluster = kmeans_find_closest_cluster(clusters, &temp_datapoint_item, data_dimension);
            if (STANDARD_ERROR_CODE == closest_cluster) {
                goto l_cleanup;
            }

            status_code = kmeans_assign_datapoint_to_cluster(clusters, &temp_datapoint_item, closest_cluster, data_dimension);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }

        /* Update centroids. */
        for (cluster_index = 0; cluster_index < clusters->item_count; cluster_index++) {
            status_code = kmeans_update_cluster_centroid(clusters, cluster_index, data_dimension, &update_centroid_diff);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }

            /* Maintain the largest centroid diff to check convergence. */
            if (update_centroid_diff > maximum_convergence_distance) {
                maximum_convergence_distance = update_centroid_diff;
            }
        }

        current_iteration++;
    } while (maximum_convergence_distance >= convergence_threshold);

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}
