import subprocess
import os
import sys
import numpy as np
from sklearn.datasets import make_blobs

# --- Setup ---
def run_command(command, goal_name):
    """
    Executes a shell command and prints the output.
    Handles potential errors in the subprocess call.
    """
    print("-" * 50)
    print(f"Running test for '{goal_name}'...")
    print(f"Command: {' '.join(command)}")
    print("-" * 50)
    try:
        # Use subprocess.run to execute the command and capture output
        result = subprocess.run(
            command, 
            capture_output=True, 
            text=True, 
            check=True
        )
        print("Output:\n")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("An Error Has Occurred")
        print(f"Error executing command: {e}")
        print(f"Stderr: {e.stderr}")
    except FileNotFoundError:
        print("An Error Has Occurred")
        print(f"Could not find executable. Make sure your C code is compiled and the Python C extension is built.")
    print("\n\n")

def create_random_data_file(file_name, N, d):
    """Creates a large dummy data file with no inherent clusters."""
    test_data = np.random.rand(N, d)
    with open(file_name, "w") as f:
        for row in test_data:
            f.write(",".join(map(str, row)) + "\n")
    print(f"Created random data file: {file_name} with {N} data points of dimension {d}")
    return file_name

def create_clusterable_data_file(file_name, N, d, centers, random_state=42):
    """Creates a data file with clear, well-defined clusters."""
    test_data, _ = make_blobs(n_samples=N, n_features=d, centers=centers, random_state=random_state)
    with open(file_name, "w") as f:
        for row in test_data:
            f.write(",".join(map(str, row)) + "\n")
    print(f"Created clusterable data file: {file_name} with {N} data points, {d} dimensions, and {centers} centers")
    return file_name

# --- Main Test Execution ---
def main():
    # --- Test with random data (stress test) ---
    print("--- Running tests with random data (low expected silhouette score) ---")
    random_file_name = create_random_data_file("random_data.txt", N=1000, d=50)

    # Test symnmf.py entry points with random data
    k_test = 3 # Choose a reasonable k for the stress test
    run_command([sys.executable, "symnmf.py", str(k_test), "sym", random_file_name], "sym (Python, Random Data)")
    run_command([sys.executable, "symnmf.py", str(k_test), "ddg", random_file_name], "ddg (Python, Random Data)")
    run_command([sys.executable, "symnmf.py", str(k_test), "norm", random_file_name], "norm (Python, Random Data)")
    run_command([sys.executable, "symnmf.py", str(k_test), "symnmf", random_file_name], "symnmf (Python, Random Data)")

    # Test symnmf C program entry points with random data
    run_command(["./symnmf", "sym", random_file_name], "sym (C, Random Data)")
    run_command(["./symnmf", "ddg", random_file_name], "ddg (C, Random Data)")
    run_command(["./symnmf", "norm", random_file_name], "norm (C, Random Data)")

    # Test analysis.py with random data
    run_command([sys.executable, "analysis.py", str(k_test), random_file_name], "analysis.py (Random Data)")
    os.remove(random_file_name)


    # --- Test with clusterable data (high expected silhouette score) ---
    print("--- Running tests with clusterable data (high expected silhouette score) ---")
    k_cluster = 3
    cluster_file_name = create_clusterable_data_file("cluster_data.txt", N=200, d=2, centers=k_cluster)

    # Test symnmf.py entry points with clusterable data
    run_command([sys.executable, "symnmf.py", str(k_cluster), "sym", cluster_file_name], "sym (Python, Clusterable Data)")
    run_command([sys.executable, "symnmf.py", str(k_cluster), "ddg", cluster_file_name], "ddg (Python, Clusterable Data)")
    run_command([sys.executable, "symnmf.py", str(k_cluster), "norm", cluster_file_name], "norm (Python, Clusterable Data)")
    run_command([sys.executable, "symnmf.py", str(k_cluster), "symnmf", cluster_file_name], "symnmf (Python, Clusterable Data)")

    # Test symnmf C program entry points with clusterable data
    run_command(["./symnmf", "sym", cluster_file_name], "sym (C, Clusterable Data)")
    run_command(["./symnmf", "ddg", cluster_file_name], "ddg (C, Clusterable Data)")
    run_command(["./symnmf", "norm", cluster_file_name], "norm (C, Clusterable Data)")

    # Test analysis.py with clusterable data
    run_command([sys.executable, "analysis.py", str(k_cluster), cluster_file_name], "analysis.py (Clusterable Data)")
    os.remove(cluster_file_name)

    print("\nAll tests completed.")

if __name__ == "__main__":
    main()
