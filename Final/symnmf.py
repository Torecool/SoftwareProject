import symnmfmodule
import sys
import math
import numpy as np

# Set random seed globally.
np.random.seed(1234)

# Constants for SymNMF behavior and validation
MINIMUM_VALID_K_VALUE = 1

# Num required args including program name.
NUM_REQUIRED_ARGS = 4
K_VALUE_ARG_INDEX = 1
GOAL_ARG_INDEX = 2
FILENAME_ARG_INDEX = 3

ERROR_EXIT_CODE = 1

# Error messages
GENERAL_ERROR_MESSAGE = 'An Error Has Occurred'


def symnmf_goal_handler(k_value, datapoints):
    """
    Perform SymNMF and output H matrix.

    Parameters
    ----------
    k_value : int
        The number of required clusters.
    datapoints : np.ndarray
        A 1D NumPy array representing the input datapoints.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the output association matrix.
    """
    # K value is used here so need to validate.
    assert MINIMUM_VALID_K_VALUE < k_value

    n_value = len(datapoints)

    # Calculate normalized similarity W matrix.
    w_matrix = symnmfmodule.norm(datapoints)

    # Average of all values of the W matrix.
    w_avg = np.mean(w_matrix)
    
    # Initialize decomposition H matrix.
    initial_h_matrix = np.random.uniform(
        low=0,
        high=2 * math.sqrt(w_avg / k_value),
        size=(n_value, k_value)
    )

    # Perform SymNMF algorithm.
    return symnmfmodule.symnmf(initial_h_matrix, w_matrix)


def sym_goal_handler(k_value, datapoints):
    """
    Calculate and output the similarity matrix.

    Parameters
    ----------
    k_value : int
        The number of required clusters.
    datapoints : np.ndarray
        A 1D NumPy array representing the input datapoints.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the output similarity matrix.
    """
    return symnmfmodule.sym(datapoints)


def ddg_goal_handler(k_value, datapoints):
    """
    Calculate and output the diagonal degree matrix.

    Parameters
    ----------
    k_value : int
        The number of required clusters.
    datapoints : np.ndarray
        A 1D NumPy array representing the input datapoints.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the output diagonal degree matrix.
    """
    return symnmfmodule.ddg(datapoints)


def norm_goal_handler(k_value, datapoints):
    """
    Calculate and output the normalized similarity matrix.

    Parameters
    ----------
    k_value : int
        The number of required clusters.
    datapoints : np.ndarray
        A 1D NumPy array representing the input datapoints.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the output normalized similarity matrix.
    """
    return symnmfmodule.norm(datapoints)


# Goals supported by the program
PROGRAM_GOALS = {
    # Perform SymNMF and output H matrix.
    'symnmf': symnmf_goal_handler,

    # Calculate and output the similarity matrix.
    'sym': sym_goal_handler,
    
    # Calculate and output the diagonal degree matrix.
    'ddg': ddg_goal_handler,
    
    # Calculate and output the normalized similarity matrix.
    'norm': norm_goal_handler
}


def parse_args():
    """
    Parse command line arguments. Should be in the form:
    python3 symnmf.py <k_param> <goal_param> <filename_param>

    Returns
    -------
    tuple
        The parsed arguments: (K, goal, filename)
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
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    try:
        goal_param = sys.argv[GOAL_ARG_INDEX]
        assert goal_param in PROGRAM_GOALS.keys()
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    # Assume file exists and is valid.
    try:
        filename_param = sys.argv[FILENAME_ARG_INDEX]
    except:
        print(GENERAL_ERROR_MESSAGE)
        raise

    return (k_param, goal_param, filename_param)


def print_results(output_matrix):
    """
    Print the selected goal's output matrix. 

    Parameters
    ----------
    output_matrix : np.ndarray
        A 2D NumPy array representing the output matrix.
    """
    for row in output_matrix:
        print(",".join(f"{val:.4f}" for val in row))


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
    return np.loadtxt(filename, delimiter=",", ndmin=2)

    
if __name__ == "__main__":
    try:
        k_param, goal_param, filename_param = parse_args()
    except:
        sys.exit(ERROR_EXIT_CODE)

    try:
        # Parse input and execute algorithm.
        datapoints = process_input(filename_param)
        if k_param >= len(datapoints):
            print(GENERAL_ERROR_MESSAGE)
            sys.exit(ERROR_EXIT_CODE)

        goal_handler = PROGRAM_GOALS[goal_param]
        output_matrix = goal_handler(k_param, datapoints)
        print_results(output_matrix)
    except Exception as e:
        print(e)
        print(GENERAL_ERROR_MESSAGE)
        sys.exit(ERROR_EXIT_CODE)
