#ifndef __KMEANS_H__
#define __KMEANS_H__

#include <stdlib.h>

/********************** Constants ******************************/
#define MINIMUM_VALID_K_VALUE (1)
#define MINIMUM_VALID_NUM_ITERATIONS (1)
#define MAXIMUM_VALID_NUM_ITERATIONS (1000)

#define DATAPOINTS_VECTOR_INITIAL_ITEM_RESERVE (32)

#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)

/********************** Structures ****************************/
/* Single datapoint item to be stored in the datapoints vector. */
struct datapoint_item {
    double *coords;
};


/* Single cluster item to be stored in the clusters vector. */
struct cluster_item {
    double *centroid;
    struct vector *datapoints;
};

/********************** Functions *****************************/
/**
 * Release the clusters vector and all its associated allocations.
 * 
 * @param clusters The clusters vector to free.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int KMEANS_free_cluster_vector(struct vector *clusters);

/**
 * Release the datapoints vector and all its associated allocations.
 * 
 * @param datapoints The datapoints vector to free.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int KMEANS_free_datapoint_vector(struct vector *datapoints);

/**
 * Initialize a cluster item within a given cluster vector.
 * 
 * @param clusters The clusters vector.
 * @param initial_centroid The initial centroid of the cluster.
 * @param data_dimension The dimension of the datapoints.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int KMEANS_create_cluster_item(
    struct vector *clusters, 
    struct datapoint_item *initial_centroid,
    size_t data_dimension
);

/**
 * Execute the K-means algorithm with the given configurations.
 * 
 * @param clusters The empty clusters vector.
 * @param datapoints The preloaded datapoints vector.
 * @param k_value The K parameter for the algorithm (cluster count).
 * @param num_iterations The iter parameter for the algorithm (maximum iterations).
 * @param data_dimension The dimension of the datapoints.
 * @param convergence_threshold The threshold value to stop iterations.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int KMEANS_perform_kmeans_clustering(
    struct vector *clusters,
    struct vector *datapoints,
    size_t k_value,
    size_t num_iterations,
    size_t data_dimension,
    double convergence_threshold
);

#endif /* __KMEANS_H__ */
