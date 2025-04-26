#!/bin/bash

# TO RUN THIS SCRIPT:
# 1. Make sure you have the C implementation (kmeans.c) and Python implementation (kmeans.py) in the same directory as this script.
# 2. Make sure you have the test input+output files in a directory named "tests" in the same directory as this script.
# 3. Run "chmod +x test.sh" to make this script executable (only for the first time, not needed later).
# 4. Run "./test.sh".

echo "Compiling C implementation..."
gcc -ansi -Wall -Wextra -Werror -pedantic-errors ./kmeans.c -o kmeans -lm
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi
echo "Compilation successful."
echo ""

# Create output directory
mkdir -p test_output

run_test() {
    test_num=$1
    k=$2
    max_iter=$3

    echo "Running test ${test_num} (K=${k}, max_iter=${max_iter})..."

    # Run C implementation
    echo "| Running C implementation..."
    ./kmeans $k $max_iter < tests/input_${test_num}.txt > test_output/c_output_${test_num}.txt
    c_status=$?

    # Run Python implementation
    # echo "| Running Python implementation..."
    # python3 ./kmeans.py $k $max_iter < tests/input_${test_num}.txt > test_output/py_output_${test_num}.txt
    # py_status=$?

    # Check for execution errors
    if [ $c_status -ne 0 ]; then
        echo ""
        echo "| C implementation failed with exit code ${c_status}"
        return 1
    fi

    # if [ $py_status -ne 0 ]; then
    #     echo ""
    #     echo "| Python implementation failed with exit code ${py_status}"
    #     return 1
    # fi

    # Compare outputs
    # diff test_output/c_output_${test_num}.txt test_output/py_output_${test_num}.txt > /dev/null
    # if [ $? -ne 0 ]; then
    #     echo ""
    #     echo "FAIL: C and Python outputs differ for test ${test_num}"
    #     return 1
    # fi

    diff test_output/c_output_${test_num}.txt tests/output_${test_num}.txt > /dev/null
    if [ $? -ne 0 ]; then
        echo ""
        echo "FAIL: Output doesn't match expected output for test ${test_num}"
        return 1
    fi

    echo "PASS: Test ${test_num} successful"
    echo ""
    return 0
}

run_test 1 3 600
test1_result=$?

run_test 2 7 ""
test2_result=$?

run_test 3 15 300
test3_result=$?


if [ $test1_result -eq 0 ] && [ $test2_result -eq 0 ] && [ $test3_result -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
