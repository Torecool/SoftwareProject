import mykmeanssp
import sys
import math
import numpy as np
import pandas as pd


# Set random seed globally.
np.random.seed(1234)

# Constants for K-means behavior and validation
DEFAULT_NUM_ITERATIONS = 400
CONVERGENCE_THRESHOLD = 0.001
MINIMUM_VALID_K_VALUE = 1
MINIMUM_VALID_NUM_ITERATIONS = 1
MAXIMUM_VALID_NUM_ITERATIONS = 1000
MINIMUM_VALID_EPSILON_VALUE = 0

# Num required args including program name.
NUM_REQUIRED_ARGS = 5
NUM_OPTIONAL_ARGS = 1
K_VALUE_ARG_INDEX = 1
NUM_ITERATIONS_ARG_INDEX = 2
EPSILON_ARG_INDEX = 3
FILENAME_1_ARG_INDEX = 4
FILENAME_2_ARG_INDEX = 5

ERROR_EXIT_CODE = 1

# Error messages
INVALID_K_VALUE_MESSAGE = 'Invalid number of clusters!'
INVALID_NUM_ITERATIONS_MESSAGE = 'Invalid maximum iteration!'
INVALID_EPSILON_MESSAGE = 'Invalid epsilon!'
GENERAL_ERROR_MESSAGE = 'An Error Has Occurred'


def find_nearest_centroid(datapoint, centeroids):
    distances = [np.linalg.norm(datapoint[1] - centroid[1]) for centroid in centeroids]
    return min(distances)


def choose_centroid_candidate(centeroid_candidates):
    distances_sum = sum(candidate[1] for candidate in centeroid_candidates)
    choices = np.array([candidate[0] for candidate in centeroid_candidates], dtype=object)
    probabilities = [candidate[1] / distances_sum for candidate in centeroid_candidates]

    return choices[np.random.choice(len(choices), p=probabilities)]


def print_results(initial_centroids, final_centroids):
    print(",".join(str(int(centroid_key)) for centroid_key, centroid_coord in initial_centroids))

    for centroid in final_centroids:
        print(','.join(f"{coord:.4f}" for coord in centroid))

        
# Runs the k-means clustering algorithm
def kmeans(datapoints, k_param, iter_param, epsilon_param):
    # Choose an initial centroid at random.
    initial_centroid = list(datapoints.items())[np.random.choice(len(datapoints))]
    del datapoints[initial_centroid[0]]

    selected_centroids = [initial_centroid]
    
    # Select remaining K - 1 centeroids.
    for iteration in range(k_param - 1):
        centeroid_candidates = []
        
        for datapoint in datapoints.items():
            # For each remaining datapoint, find the nearest centeroid and the distance to it.
            centeroid_candidates.append((datapoint, find_nearest_centroid(datapoint, selected_centroids)))
            
        # Choose the new centeroid according to weighted probability distribution.
        chosen_candidate = choose_centroid_candidate(centeroid_candidates)
        selected_centroids.append(chosen_candidate)
        del datapoints[chosen_candidate[0]]

    centroids_coords = [tuple(centroid[1]) for centroid in selected_centroids]
    datapoints_coords = [tuple(datapoint) for datapoint in datapoints.values()]
    return selected_centroids, mykmeanssp.fit(centroids_coords, datapoints_coords, iter_param, epsilon_param)


def parse_args():
    if len(sys.argv) < NUM_REQUIRED_ARGS or len(sys.argv) > NUM_REQUIRED_ARGS + NUM_OPTIONAL_ARGS:
        print(GENERAL_ERROR_MESSAGE)
        raise

    try:
        # Validate mandatory parameters
        k_param = float(sys.argv[K_VALUE_ARG_INDEX])
        if k_param.is_integer():
            k_param = int(k_param)
        else:
            raise ValueError(f"Expected an integer-like value, got non-integer: {k_param}")
        
        assert MINIMUM_VALID_K_VALUE < k_param
    except:
        print(INVALID_K_VALUE_MESSAGE)
        raise

    optional_parameter_shift = 0
    try:
        # Validate optional parameters
        if len(sys.argv) < NUM_REQUIRED_ARGS + NUM_OPTIONAL_ARGS:
            iter_param = DEFAULT_NUM_ITERATIONS
            optional_parameter_shift = 1
        else:
            iter_param = float(sys.argv[NUM_ITERATIONS_ARG_INDEX])
            if iter_param.is_integer():
                iter_param = int(iter_param)
            else:
                raise ValueError(f"Expected an integer-like value, got non-integer: {iter_param}")
                
        assert MINIMUM_VALID_NUM_ITERATIONS < iter_param and iter_param < MAXIMUM_VALID_NUM_ITERATIONS
    except:
        print(INVALID_NUM_ITERATIONS_MESSAGE)
        raise

    # Validate more mandatory parameters
    try:
        epsilon_param = float(sys.argv[EPSILON_ARG_INDEX - optional_parameter_shift])
        assert MINIMUM_VALID_EPSILON_VALUE <= epsilon_param
    except:
        print(INVALID_EPSILON_MESSAGE)
        raise

    try:
        filename_1_param = sys.argv[FILENAME_1_ARG_INDEX - optional_parameter_shift]
        filename_2_param = sys.argv[FILENAME_2_ARG_INDEX - optional_parameter_shift]
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    return (k_param, iter_param, epsilon_param, filename_1_param, filename_2_param)


def process_input(filename_1, filename_2):
    # Read the CSV files
    df1 = pd.read_csv(filename_1, header=None)
    df2 = pd.read_csv(filename_2, header=None)

    # Rename first column to 'key' for both
    df1.rename(columns={0: 'key'}, inplace=True)
    df2.rename(columns={0: 'key'}, inplace=True)

    # Perform inner join on 'key' column
    joined_df = pd.merge(df1, df2, on='key', how='inner')

    # Sort in ascending order.
    joined_df = joined_df.sort_values(by='key')

    # Extract datapoints.
    datapoints = {
        row[0]: np.array(row[1:])
        for row in joined_df.to_numpy()
    }

    return datapoints

    
if __name__ == "__main__":
    try:
        k_param, iter_param, epsilon_param, filename_1_param, filename_2_param = parse_args()
    except:
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Parse input and execute algorithm.
        datapoints = process_input(filename_1_param, filename_2_param)
        if k_param >= len(datapoints):
            print(INVALID_K_VALUE_MESSAGE)
            sys.exit(ERROR_EXIT_CODE)

        selected_centroids, final_centroids = kmeans(datapoints, k_param, iter_param, epsilon_param)
        print_results(selected_centroids, final_centroids)
    except Exception as e:
        print(GENERAL_ERROR_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)
