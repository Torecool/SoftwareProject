import subprocess
import os
import sys
import numpy as np

# --- Setup ---
# Create a dummy data file for testing.
# The project document mentions a .txt file with N data points.
# We'll create a simple 5x3 matrix for our test.
test_data = [
    [1.0, 2.0, 3.0],
    [1.1, 2.1, 3.1],
    [5.0, 6.0, 7.0],
    [5.1, 6.1, 7.1],
    [1.5, 2.5, 3.5]
]
file_name = "test_data.txt"

with open(file_name, "w") as f:
    for row in test_data:
        f.write(",".join(map(str, row)) + "\n")

print(f"Created dummy data file: {file_name}\n")


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

# --- symnmf.py tests ---
# Test the 'sym' goal
k_sym = 2
goal_sym = "sym"
run_command([sys.executable, "symnmf.py", str(k_sym), goal_sym, file_name], "sym (Python)")

# Test the 'ddg' goal
k_ddg = 2
goal_ddg = "ddg"
run_command([sys.executable, "symnmf.py", str(k_ddg), goal_ddg, file_name], "ddg (Python)")

# Test the 'norm' goal
k_norm = 2
goal_norm = "norm"
run_command([sys.executable, "symnmf.py", str(k_norm), goal_norm, file_name], "norm (Python)")

# Test the 'symnmf' goal
k_nmf = 2
goal_nmf = "symnmf"
run_command([sys.executable, "symnmf.py", str(k_nmf), goal_nmf, file_name], "symnmf (Python)")

# --- symnmf (C program) tests ---
# Test the 'sym' goal for the C program
run_command(["./symnmf", "sym", file_name], "sym (C)")

# Test the 'ddg' goal for the C program
run_command(["./symnmf", "ddg", file_name], "ddg (C)")

# Test the 'norm' goal for the C program
run_command(["./symnmf", "norm", file_name], "norm (C)")

# --- analysis.py test ---
# Test the analysis script
run_command([sys.executable, "analysis.py", str(k_nmf), file_name], "analysis.py")

# --- Cleanup ---
# Remove the temporary file after the tests are done
os.remove(file_name)
print(f"Removed dummy data file: {file_name}")
print("All tests completed.")
