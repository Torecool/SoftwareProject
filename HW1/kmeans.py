import sys
import math

# Constants for K-means behavior and validation
DEFAULT_NUM_ITERATIONS = 400
CONVERGENCE_THRESHOLD = 0.001
MINIMUM_VALID_K_VALUE = 1
MINIMUM_VALID_NUM_ITERATIONS = 1
MAXIMUM_VALID_NUM_ITERATIONS = 1000

# Magic value representing that no cluster is close to the given datapoint.
MAX_DIST = 1000

# Num required args including program name.
NUM_REQUIRED_ARGS = 2
NUM_OPTIONAL_ARGS = 1
K_VALUE_ARG_INDEX = 1
NUM_ITERATIONS_ARG_INDEX = 2

ERROR_EXIT_CODE = 1

# Error messages
INVALID_K_VALUE_MESSAGE = 'Incorrect number of clusters!'
INVALID_NUM_ITERATIONS_MESSAGE = 'Incorrect maximum iteration!'
GENERAL_ERROR_MESSAGE = 'An Error Has Occurred'


# Cluster class represents a single cluster in the k-means algorithm
class Cluster:
    def __init__(self, centroid, datapoints):
        self.centroid = centroid          # The current centroid (list of floats)
        self.datapoints = datapoints      # List of data points assigned to this cluster

    def __str__(self):
        # Returns a string representation of the centroid, formatted to 4 decimal places
        return ','.join(f"{coord:.4f}" for coord in self.centroid)

    def get_centroid(self):
        # Returns the centroid of the cluster
        return self.centroid

    def get_datapoints(self):
        # Returns the list of datapoints assigned to the cluster
        return self.datapoints

    def get_size(self):
        # Returns the number of datapoints in the cluster
        return len(self.datapoints)

    def add_datapoint(self, dp):
        # Adds a datapoint to the cluster
        self.datapoints.append(dp)

    def empty(self):
        # Removes all datapoints from the cluster (used before reassigning in next iteration)
        self.datapoints = []

    def update_centroid(self):
        # Recalculates the centroid based on the current datapoints
        if len(self.datapoints) > 0:
            for i in range(len(self.centroid)):
                # Update each coordinate of the centroid as the mean of that coordinate over all datapoints
                self.centroid[i] = sum(vec[i] for vec in self.datapoints) / len(self.datapoints)

    def calc_dist(self, dp):
        # Calculates the Euclidean distance from a datapoint to the cluster's centroid
        temp_sum = 0
        for i in range(len(self.centroid)):
            temp_sum += (self.centroid[i] - dp[i]) ** 2
        return math.sqrt(temp_sum)


# Reads input data from standard input (stdin) and returns a list of datapoints
def parse_input():
    data = sys.stdin.read().strip().split('\n')
    data_vec = [tuple(map(float, line.strip().split(','))) for line in data if line]
    return data_vec


# Assigns a datapoint to the cluster with the nearest centroid
def assign_datapoint(datapoint, clusters_list):
    min_distance = MAX_DIST
    min_cluster = None

    for cluster in clusters_list:
        dist = cluster.calc_dist(datapoint)
        if dist < min_distance:
            min_distance = dist
            min_cluster = cluster

    min_cluster.add_datapoint(datapoint)
    return


# Initializes k clusters by selecting the first k datapoints as initial centroids
def initialize_kmeans(data_vec, k_param):
    if len(data_vec) <= k_param:
        raise ValueError()

    clusters_list = []
    for i in range(k_param):
        # Initialize cluster with centroid
        clusters_list.append(Cluster(centroid=list(data_vec[i]), datapoints=[]))

    return clusters_list


# Runs the k-means clustering algorithm
def kmeans(data_vec, k_param, iter_param=DEFAULT_NUM_ITERATIONS):
    # Initialize data and clusters
    clusters_list = initialize_kmeans(data_vec, k_param)

    for curr_iter in range(iter_param):
        
        # Clear all clusters for the next iteration
        for cluster in clusters_list:
            cluster.empty()

        # Save current centroids to check for convergence later
        prev_centroids = [cluster.get_centroid().copy() for cluster in clusters_list]

        # Assign all datapoints to the nearest cluster
        for datapoint in data_vec:
            assign_datapoint(datapoint, clusters_list)

        # Update centroids based on assigned datapoints
        for cluster in clusters_list:
            cluster.update_centroid()

        # Check for convergence
        cluster_close_enough = True
        for cluster, prev_centroid in zip(clusters_list, prev_centroids):
            if cluster.calc_dist(prev_centroid) >= CONVERGENCE_THRESHOLD:
                cluster_close_enough = False
                # Stop checking convergence if any centroid moved significantly
                break 

        if cluster_close_enough:
            break

    for cluster in clusters_list:
        print(cluster)


if __name__ == "__main__":
    if len(sys.argv) < NUM_REQUIRED_ARGS:
        print(GENERAL_ERROR_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Validate mandatory parameters
        k_param = int(sys.argv[K_VALUE_ARG_INDEX])
        assert MINIMUM_VALID_K_VALUE < k_param and k_param < MAXIMUM_VALID_NUM_ITERATIONS
    except:
        print(INVALID_K_VALUE_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Validate optional parameters
        if len(sys.argv) < NUM_REQUIRED_ARGS + NUM_OPTIONAL_ARGS:
            iter_param = DEFAULT_NUM_ITERATIONS
        else:
            iter_param = int(sys.argv[NUM_ITERATIONS_ARG_INDEX])

        assert MINIMUM_VALID_NUM_ITERATIONS < iter_param and iter_param < MAXIMUM_VALID_NUM_ITERATIONS
    except:
        print(INVALID_NUM_ITERATIONS_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Parse input and execute algorithm.
        data_vec = parse_input()
        kmeans(data_vec, k_param, iter_param)
    except:
        print(GENERAL_ERROR_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)
