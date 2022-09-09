import sys
from tabnanny import verbose

import parse
import problem as pb
import checker
from solver import *


def mainFunction(
    instanceFilename: str, assignmentFilename: str, verbose: bool, maxTime: int
):
    if verbose:
        print("Received instance files ", instanceFilename, " and ", assignmentFilename)

    # reading the data
    data = parse.parseFiles(instanceFilename, assignmentFilename)

    if verbose:
        pb.printData(data)

    solutionTemp = pb.Solution(data.initialAssignment, 0)
    checker.check(data, solutionTemp, verbose)  # check if the solution is feasible

    # solve the problem
    solution = solve(data, maxTime, verbose=False)

    if verbose:  # print the solution
        print("Solution of cost ", solution.cost)
        for i in range(data.nbProcess):
            print(
                solution.assignment[i], end=" "
            )  # here the solution is written with our convention 1..n
        print()

    # test if the solution is feasible and compute its cost
    checker.check(data, solution, verbose)  # check if the solution is feasible


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "You should provide exactly 4 arguments: instance_filename, assignment_filename, time_limit (s) and verbose (int)"
        )
        print(
            "Usage: python3 src/run.py <model_file> <assignment_file> <time_limit (s)> <verbose (int)>"
        )
        exit(1)

    instance_filename = sys.argv[1]
    assignment_filename = sys.argv[2]
    time_limit = int(sys.argv[3])
    verbose = bool(int(sys.argv[4]))

    mainFunction(instance_filename, assignment_filename, verbose, time_limit)
