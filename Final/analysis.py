from symnmf import symnmf_goal_handler
from kmeans import kmeans
from sklearn.metrics import silhouette_score
import numpy as np
import sys

# Constants for behavior and validation
MINIMUM_VALID_K_VALUE = 1
MAX_ITERATIONS = 300

# Num required args including program name.
NUM_REQUIRED_ARGS = 3
K_VALUE_ARG_INDEX = 1
FILENAME_ARG_INDEX = 2

ERROR_EXIT_CODE = 1

# Error messages
GENERAL_ERROR_MESSAGE = 'An Error Has Occurred'


def parse_args():
    """
    Parse command line arguments. Should be in the form:
    python3 analysis.py <k_param> <filename_param>

    Returns
    -------
    tuple
        The parsed arguments: (K, filename)
    """
    if len(sys.argv) != NUM_REQUIRED_ARGS:
        print(GENERAL_ERROR_MESSAGE)
        raise ValueError(f"Invalid number of arguments: {len(sys.argv)}")

    # Validate mandatory parameters
    try:
        k_param = float(sys.argv[K_VALUE_ARG_INDEX])
        if k_param.is_integer():
            k_param = int(k_param)
        else:
            raise ValueError(f"Expected an integer-like value, got non-integer: {k_param}")
        
        assert MINIMUM_VALID_K_VALUE < k_param
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    # Assume file exists and is valid.
    try:
        filename_param = sys.argv[FILENAME_ARG_INDEX]
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    return (k_param, filename_param)


def print_results(score_kmeans, score_symnmf):
    """
    Print out the results of the analysis.

    Parameters
    ----------
    score_kmeans : float
        The silhouette coefficient of the K-Means algorithm.
    score_symnmf : float
        The silhouette coefficient of the SymNMF algorithm.
    """
    print(f"sym: {score_symnmf}")
    print(f"kmeans: {score_kmeans}")


def process_input(filename):
    """
    Read datapoints from an input file. 

    Parameters
    ----------
    filename : str
        The filename containing the datapoints.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the input datapoints.
    """
    return np.loadtxt(filename, delimiter=",")


def perform_analysis(k_value, datapoints):
    # Execute kmeans algorithm.
    datapoint_list = [tuple(point) for point in datapoints]
    clusters_list = kmeans(datapoint_list, k_value, MAX_ITERATIONS)

    # Extract labels from cluster output.
    kmeans_labels = np.empty(len(datapoints), dtype=int)
    for datapoint_index, datapoint in enumerate(datapoints):
        for cluster_index, cluster in enumerate(clusters_list):
            if tuple(datapoint) in cluster.datapoints:
                kmeans_labels[datapoint_index] = cluster_index

    # Execute SymNMF algorithm.
    final_h_matrix = symnmf_goal_handler(k_value, datapoints)
    # Assign labels by taking argmax of rows in final H.
    symnmf_labels = np.argmax(final_h_matrix, axis=1)

    # Compute silhouette scores.
    score_kmeans = silhouette_score(datapoints, kmeans_labels)
    score_symnmf = silhouette_score(datapoints, symnmf_labels)

    return (score_kmeans, score_symnmf)


if __name__ == "__main__":
    try:
        k_param, filename_param = parse_args()
    except:
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Parse input perform analysis.
        datapoints = process_input(filename_param)
        if k_param >= len(datapoints):
            print(GENERAL_ERROR_MESSAGE)
            sys.exit(ERROR_EXIT_CODE)

        (score_kmeans, score_symnmf) = perform_analysis(k_param, datapoints)
        print_results(score_kmeans, score_symnmf)
    except Exception as e:
        print(GENERAL_ERROR_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)
