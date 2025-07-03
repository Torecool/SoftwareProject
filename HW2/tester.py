#!/usr/bin/env python3
import argparse
import decimal
from decimal import Decimal
import io
import itertools
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
from typing import Callable, Optional, Sequence

import numpy as np

__author__ = "Yahel Caspi"
__version__ = "1.2.1"

LINE_CONFIG_REGEX = re.compile(
    r"(?P<idx>\d+)\. k=(?P<k>\d+), max_iter\s+=\s+(?:not provided|(?P<max_iter>\d+)),\s+eps=(?P<eps>\d+(?:\.\d+)?),\s+(?P<filename1>\w+),\s+(?P<filename2>\w+)"
)
VALGRIND_ERRCODE = 99

IGNORE_ERRORCODE_0 = False
USE_VALGRIND = False


def print_green(msg: str):
    print(f"\033[32m{msg}\033[0m")


def print_yellow(msg: str):
    print(f"\033[33m{msg}\033[0m")


def print_red(msg: str):
    print(f"\033[31m{msg}\033[0m")


def print_white_on_red(msg: str):
    print(f"\033[97;41m{msg}\033[0m")


def generate_data(K=20, points_num=None):
    rng = np.random.default_rng()
    dim = rng.integers(2, 10)
    N = points_num if points_num else rng.integers(100, 700)

    centroids = rng.uniform(-11, 11, (K, dim))
    data = rng.choice(centroids, N) + rng.standard_normal((N, dim))

    # Remove duplicate points
    data = np.unique(data, axis=0)

    while data.shape[0] < N:
        extra = rng.choice(centroids, N - data.shape[0]) + rng.standard_normal(
            (N - data.shape[0], dim)
        )
        data = np.unique(np.vstack([data, extra]), axis=0)

    return data.astype(np.float64)


def split_data(data: np.ndarray):
    rng = np.random.default_rng()
    ids = np.arange(data.shape[0], dtype=np.float64).reshape(-1, 1)

    split_at = rng.integers(1, data.shape[1] - 1) if data.shape[1] > 2 else 1
    part1, part2 = np.hsplit(data, [split_at])
    part1 = np.hstack((ids, part1))
    part2 = np.hstack((ids, part2))

    rng.shuffle(part1)
    rng.shuffle(part2)

    return part1, part2


def make_stub_files():
    file1 = tempfile.NamedTemporaryFile()
    file2 = tempfile.NamedTemporaryFile()
    buf1 = io.TextIOWrapper(file1)
    buf2 = io.TextIOWrapper(file2)

    data = generate_data(points_num=100)
    data1, data2 = split_data(data)

    for buf, data in ((buf1, data1), (buf2, data2)):
        for row in data:
            print(",".join(f"{x:.4f}" for x in row), file=buf)
        buf.detach()

    return file1, file2


def generate_invalid_param(test_eps: bool):
    english = ("bug", "a")
    hebrew = ("באג", "א")
    cjk = ("虫", "バグ")
    emojis = ("🐞", "🏴‍☠️")
    control_chars = ("\a", "\b")

    negative_numbers = ("-1", "-1.00", "-1.5", "-0.3", "-5")
    zeros = ("0", "0.0")
    positive_floats = ("1.23", "0.42")
    big_integers = ("65,539", "65539")
    ones = ("1", "1.0000")

    good_integers = ("3", "3.000", "03")

    yield from itertools.chain(
        english, hebrew, cjk, emojis, control_chars, negative_numbers
    )

    if not test_eps:
        yield from itertools.chain(big_integers, ones, positive_floats, zeros)

    yield from map(
        lambda x: "".join(x),
        itertools.product(good_integers, itertools.chain(english, hebrew, cjk, emojis)),
    )


def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True, dest="command")

    main_parser = subparsers.add_parser(
        "main", help="test the python app (with the C module)"
    )
    main_parser.add_argument(
        "-e0",
        "--ignore-errorcode-0",
        action="store_true",
        help="allows the code to exit with status 0 when an error is expected",
    )
    main_parser.add_argument(
        "--valgrind",
        action="store_true",
        help="use valgrind to test the C module for memory leaks",
    )
    main_parser.add_argument(
        "tests_dir",
        help="path to the tests directory (from Moodle)",
        metavar="TESTS_DIR",
    )

    c_parser = subparsers.add_parser(
        "c", help="test the C code - specifically, " "the fit function"
    )
    c_parser.add_argument(
        "--trials", type=int, default=10, help="set the number of trials"
    )

    tar_parser = subparsers.add_parser("tar-gz", help="verify the tar.gz file")
    tar_parser.add_argument("--id1", required=True)
    tar_parser.add_argument("--id2")
    tar_parser.add_argument(
        "--dir",
        help="the directory in which the tar.gz file is located "
        "(defaults to the current directory)",
    )

    return parser


def run_test_files(tests_dir: Path):
    # Read configuration
    readme_path = tests_dir / Path("test_readme.txt")
    configs = []
    with readme_path.open() as f:
        for line in f:
            match = LINE_CONFIG_REGEX.match(line)
            if match:
                config = match.groupdict()
                config["filename1"] = str(
                    tests_dir / Path(config["filename1"]).with_suffix(".txt")
                )
                config["filename2"] = str(
                    tests_dir / Path(config["filename2"]).with_suffix(".txt")
                )
                configs.append(config)

    for config in configs:
        success = True

        print("Test", config["idx"])
        result, valgrind_logfile = execute(config)

        if result.returncode != 0 and (
            not USE_VALGRIND or result.returncode != VALGRIND_ERRCODE
        ):
            success = False
            print(f"process returned with code {result.returncode}")
        elif USE_VALGRIND and result.returncode == VALGRIND_ERRCODE:
            valgrind_log = valgrind_logfile.read()  # type: ignore
            print_red("memory leak detected by valgrind")
            print_white_on_red(valgrind_log.decode())

        if result.stderr:
            success = False
            print("process had non-empty stderr:")
            print_white_on_red(result.stderr)

        # Compare outputs
        output_path = tests_dir / Path(f"output_{config['idx']}.txt")
        if (
            reference_output := output_path.read_text().rstrip()
        ) != result.stdout.rstrip():
            result_lines = result.stdout.rstrip().splitlines()
            reference_output_lines = reference_output.splitlines()

            if not verify_outputs(result_lines, reference_output_lines):
                success = False
                print("mismatch between the process output and the target output")

        if success:
            print_green("success")
        else:
            print_red("failure")

        if valgrind_logfile:
            valgrind_logfile.close()


def verify_outputs(result: Sequence[str], reference: Sequence[str]):
    PERMITTED_DELTA = Decimal("0.0001")

    if len(result) != len(reference):
        return False

    if result[0] != reference[0]:
        return False

    for i in range(1, len(result)):
        try:
            point = tuple(map(Decimal, result[i].split(",")))
        except decimal.InvalidOperation:
            return False

        reference_point = tuple(map(Decimal, reference[i].split(",")))

        if len(point) != len(reference_point):
            return False

        if any(
            (a - b).copy_abs() > PERMITTED_DELTA for a, b in zip(point, reference_point)
        ):
            return False

    return True


def test_input_handling():
    file1, file2 = make_stub_files()

    # Test invalid parameters for each of "k", "max_iter", and "eps"
    print("Test invalid parameters")
    params = {
        "k": {
            "valid_value": "3",
            "error_msg": "Invalid number of clusters!",
        },
        "max_iter": {
            "valid_value": "300",
            "error_msg": "Invalid maximum iteration!",
        },
        "eps": {
            "valid_value": "0.001",
            "error_msg": "Invalid epsilon!",
        },
    }

    valid_values = {param: val["valid_value"] for param, val in params.items()}

    for param_name in params:
        total_tests = 0
        passed_tests = 0

        for param_value in generate_invalid_param(test_eps=(param_name == "eps")):
            # Print in magenta, erase later
            print(f"\033[35mtesting {param_name}={param_value!r}...\033[0m", end="\r")
            total_tests += 1
            config = {
                **valid_values,
                "filename1": file1.name,
                "filename2": file2.name,
            }
            config[param_name] = param_value
            result, valgrind_logfile = execute(config)

            # Erase message
            print("\033[0K", end="")

            expected_result = True
            if result.returncode != 1 and not (
                result.returncode == 0 and IGNORE_ERRORCODE_0
                or result.returncode == VALGRIND_ERRCODE and USE_VALGRIND
            ):
                expected_result = False
                print_red(
                    f"failure: process accepted invalid {param_name} {param_value!r}"
                )
            elif USE_VALGRIND and result.returncode == VALGRIND_ERRCODE:
                valgrind_log = valgrind_logfile.read()  # type: ignore
                print_red("memory leak detected by valgrind")
                print_white_on_red(valgrind_log.decode())

            if result.stdout.rstrip("\n") != params[param_name]["error_msg"]:
                expected_result = False
                print_red(
                    f"failure: process returned an incorrect error message for invalid {param_name}"
                )

            if result.stderr:
                expected_result = False
                print_red(f"failure: process had a non-empty stderr")
                print_white_on_red(result.stderr)

            if expected_result:
                passed_tests += 1
            else:
                print(f"An error occurred with {param_name}={param_value!r}")

            if valgrind_logfile:
                valgrind_logfile.close()

        if passed_tests > 0:
            print_green(
                f"{passed_tests}/{total_tests} successes: process rejected "
                f"invalid values of {param_name}"
            )

    # Check the use of too much parameters
    config = {
        **valid_values,
        "filename1": file1.name,
        "filename2": file2.name,
        "additional_args": ("123",),
    }
    result, valgrind_logfile = execute(config)

    expected_result = True
    if result.returncode != 1 and not (result.returncode == 0 and IGNORE_ERRORCODE_0
                                       or result.returncode == VALGRIND_ERRCODE and USE_VALGRIND):
        expected_result = False
        print_red("failure: process accepted too many arguments")
    elif USE_VALGRIND and result.returncode == VALGRIND_ERRCODE:
        valgrind_log = valgrind_logfile.read()  # type: ignore
        print_red("memory leak detected by valgrind")
        print_white_on_red(valgrind_log.decode())

    if result.stdout.rstrip("\n") != "An Error Has Occurred":
        expected_result = False
        print_red(
            "failure: process returned an incorrect error message for "
            "having too many arguments"
        )

    if result.stderr:
        expected_result = False
        print_red(f"failure: process had a non-empty stderr")
        print_white_on_red(result.stderr)

    if expected_result:
        print_green("success: process rejected an execution with too many arguments")

    if valgrind_logfile:
        valgrind_logfile.close()

    # Check handling of good inputs for k, max_iter, and eps
    good_inputs = [
        {"k": "3", "max_iter": "300", "eps": "0.001"},
        {"k": "010", "max_iter": "0300", "eps": "00.001"},
        {"k": "3.00", "max_iter": "300.00", "eps": "0.001"},
    ]
    for params in good_inputs:
        config = {
            **params,
            "filename1": file1.name,
            "filename2": file2.name,
        }

        result, valgrind_logfile = execute(config)

        if result.returncode == 0 and not result.stderr:
            print_green(f"success: process accepted good input {params}")
        else:
            if USE_VALGRIND and result.returncode == VALGRIND_ERRCODE:
                valgrind_log = valgrind_logfile.read()  # type: ignore
                print_red("memory leak detected by valgrind")
                print_white_on_red(valgrind_log)
            else:
                print_red(f"failure: process rejected good input {params}")

            if result.stderr:
                print("process had a non-empty stderr:")
                print_white_on_red(result.stderr)

        if valgrind_logfile:
            valgrind_logfile.close()


def fit_adapter(
    fit: Callable,
    datapoints: np.ndarray,
    initial_centroids: np.ndarray,
    eps: float,
    max_iter: int,
) -> np.ndarray:
    """Calls the `fit` functions with the given arguments and returns centriods as a numpy array.

    Edit this function to make sure the tests work with your code
    """
    result = fit(datapoints, np.ascontiguousarray(initial_centroids), eps, max_iter)
    return np.asarray(result)


def kmeans_reference(datapoints, init_centroids, eps=0.001, max_iter=300):
    datapoints = np.asarray(datapoints)
    K = init_centroids.shape[0]
    centroids = init_centroids.copy()
    for _ in range(max_iter):
        # Assign clusters
        # NOTE: the arrays are reshaped so that the distance calculations be correct.
        dists = np.linalg.norm(datapoints[:, None, :] - centroids[None, :, :], axis=2)
        labels = np.argmin(dists, axis=1)
        new_centroids = np.zeros_like(centroids)
        clusters_sizes = np.zeros(K, dtype=int)

        for i, x in enumerate(datapoints):
            k = labels[i]
            clusters_sizes[k] += 1
            new_centroids[k] += x

        for k in range(K):
            if clusters_sizes[k] > 0:
                new_centroids[k] /= clusters_sizes[k]
            else:
                # Keep the old centroid if cluster is empty
                new_centroids[k] = centroids[k]

        # Convergence check
        if np.all(np.linalg.norm(centroids - new_centroids, axis=1) < eps):
            break

        centroids = new_centroids

    return centroids


def test_fit(trials=5):
    print("Testing 'mykmeanspp.fit'")
    try:
        from sklearn import cluster as skl_cluster
    except ModuleNotFoundError:
        print_red("This test requires scikit-learn to be installed")
        return

    try:
        import mykmeanspp
    except ModuleNotFoundError:
        print_red("'mykmeanspp' module not found")
        return

    try:
        fit = mykmeanspp.fit
    except AttributeError:
        print_red("'mykmeanspp.fit' not found")

    for trial in range(1, trials + 1):
        print(f"trial {trial}:", end=" ")
        K = np.random.randint(5, 25)
        datapoints = generate_data(K)
        initial_centroids, _ = skl_cluster.kmeans_plusplus(datapoints, K)
        eps = np.random.uniform(0, 0.001)
        max_iter = np.random.randint(200, 400)

        result = fit_adapter(fit, datapoints, initial_centroids, eps, max_iter)

        reference_centroids = kmeans_reference(
            datapoints, initial_centroids, eps, max_iter
        )

        centroid_diffs = np.abs(result - reference_centroids)

        if not np.all(centroid_diffs < eps):
            print_red(
                "failure: `fit` kmeans algorithm doesn't match that of scikit-learn"
            )
            print(eps, centroid_diffs.max())
        else:
            print_green("success")


from typing import IO, Any


def execute(config) -> tuple[subprocess.CompletedProcess[str], Optional[IO[Any]]]:
    args = ["python3", "kmeanspp.py", config["k"]]
    if config.get("max_iter"):
        args.append(config["max_iter"])
    args += [config["eps"], config["filename1"], config["filename2"]]

    pass_fds = []
    if USE_VALGRIND:
        logfile = tempfile.TemporaryFile()
        fd = logfile.fileno()
        # Ensures the fd is inherited by child processes
        os.set_inheritable(fd, True)

        args = [
            "valgrind",
            "--leak-check=full",
            f"--log-fd={fd}",
            f"--error-exitcode={VALGRIND_ERRCODE}",
            "--suppressions=python.supp",
            "--show-leak-kinds=definite,indirect",
            "--errors-for-leak-kinds=definite,indirect",
        ] + args

        pass_fds.append(fd)

    if config.get("additional_args"):
        args.extend(config["additional_args"])

    result = subprocess.run(
        args,
        capture_output=True,
        text=True,
        pass_fds=pass_fds,
    )

    if USE_VALGRIND:
        logfile.seek(0)
        return result, logfile

    return result, None


def test_tar(id1: str, id2: str = "", basepath: Optional[Path] = None):
    if not id2:
        id2 = "111111111"

    basic_name = f"{id1}_{id2}_assignment2"
    filename = Path(f"{basic_name}.tar.gz")
    filepath = filename if not basepath else basepath / filename

    if not filepath.exists():
        print_red("didn't found a matching tar.gz file")

    has_dir = False
    has_files = {
        "kmeanspp.py": False,
        "kmeansmodule.c": False,
        "setup.py": False,
    }
    has_bonus = False

    with tarfile.open(filepath) as tar:
        for tarinfo in tar.getmembers():
            unknown_file = True

            if tarinfo.isdir() and tarinfo.name == basic_name:
                has_dir = True
                unknown_file = False
                print_green(f"✅ {basic_name} directory is found")
            elif tarinfo.isfile():
                if tarinfo.name == f"{basic_name}/bonus.py":
                    has_bonus = True
                    unknown_file = False
                    print_green(f"✅ bonus.py is found ✨")
                else:
                    for filename in has_files:
                        if tarinfo.name == f"{basic_name}/{filename}":
                            has_files[filename] = True
                            unknown_file = False
                            print_green(f"✅ {filename} is found")
                            break
                    else:
                        if match := re.match(
                            rf"{basic_name}/(\w+\.[ch])$", tarinfo.name
                        ):
                            unknown_file = False
                            print_yellow(f"💡 {match.group(1)} is found")

            if unknown_file:
                print_red(f"🚨 unknown file: {tarinfo.name}")

    if not has_dir:
        print_red(f"❌ missing {basic_name} directory")
    else:
        for filename in filter(lambda x: not has_files[x], has_files):
            print_red(f"❌ missing {filename}")

        if not has_bonus:
            print("bonus.py wasn't found (but isn't required)")


def main():
    args = setup_argparser().parse_args()

    match args.command:
        case "main":
            global USE_VALGRIND, IGNORE_ERRORCODE_0
            USE_VALGRIND = args.valgrind
            IGNORE_ERRORCODE_0 = args.ignore_errorcode_0

            run_test_files(Path(args.tests_dir))
            test_input_handling()
        case "c":
            test_fit(args.trials)
        case "tar-gz":
            tgz_dir = Path(args.dir) if args.dir else None
            test_tar(args.id1, args.id2, tgz_dir)


if __name__ == "__main__":
    main()
