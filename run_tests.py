# run_tests.py
# Checks that every parallel run gives exactly the same CSV as the sequential run

import filecmp
import os
import subprocess
import sys

# the program to test (works on Windows, Linux and Mac)
if len(sys.argv) > 1:
    exe = sys.argv[1]
elif os.name == "nt":
    exe = "sim.exe"
else:
    exe = "./sim"

N = 200
days = 300
seed = 42

os.makedirs("out/tests", exist_ok=True)


def run(mode, threads, chunks, part, outfile):
    """Run the simulation once and write its CSV to outfile."""
    command = [exe,
               "--N", str(N),
               "--days", str(days),
               "--seed", str(seed),
               "--mode", mode,
               "--threads", str(threads),
               "--chunks", str(chunks),
               "--partition", part,
               "--out", outfile]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


# 1. The sequential reference
reference = "out/tests/ref.csv"
run("sequential", 1, 1, "strips", reference)

passed = 0
failed = 0

# 2. Every thread count and every way of dividing the grid
for T in [1, 2, 3, 4, 7, 8, 16]:
    for part in ["strips", "blocks", "interleaved"]:

        # each case is: mode, threads, chunks
        cases = [
            ("thread-per-day", T, T),
            ("pool", T, T),
            ("pool", T, 4 * T),
        ]

        for mode, threads, chunks in cases:
            outfile = f"out/tests/{mode}_T{threads}_C{chunks}_{part}.csv"
            run(mode, threads, chunks, part, outfile)

            # compare the file contents byte by byte
            if filecmp.cmp(reference, outfile, shallow=False):
                passed += 1
            else:
                failed += 1
                print("FAIL:", outfile)

print(f"Passed: {passed}   Failed: {failed}")