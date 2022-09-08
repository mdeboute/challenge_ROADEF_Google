import sys
from tabnanny import verbose

import parse
import problem as pb
import checker
from solver import *


def mainFunction(
    instanceFilename: str, assignmentFilename: str, verbose: bool, maxTime: int
):
    print("Received instance files ", instanceFilename, " and ", assignmentFilename)

    # reading the data
    data = parse.parseFiles(instanceFilename, assignmentFilename)

    if verbose:
        pb.printData(data)

    solutionTemp = pb.Solution(data.initialAssignment, 0)
    checker.check(data, solutionTemp, verbose)  # check if the solution is feasible

    # solve the problem
    solution = solve(data, maxTime, verbose)

    if verbose:  # print the solution
        print("Solution of cost ", solution.cost)
        for i in range(data.nbProcess):
            print(
                solution.assignment[i], end=" "
            )  # here the solution is written with our convention 1..n
        print()

    # test if the solution is feasible and compute its cost
    checker.check(data, solution, True)  # check if the solution is feasible


if __name__ == "__main__":
    instanceFilename = "./data/dataA/model_a1_3.txt"
    assignmentFilename = "./data/dataA/assignment_a1_3.txt"
    timeLimit = 600
    verbose = True

    if len(sys.argv) == 4:
        print(
            "You should provide exactly 4 arguments : instanceFilename, assignmentFilename, timeLimit and verbose"
        )
        print(
            "Usage: python3 run.py data/dataA/model_a1_1.txt data/dataA/assignment_a1_1.txt 300 True"
        )

        instanceFilename = sys.argv[1]
        assignmentFilename = sys.argv[2]
        timeLimit = int(sys.argv[3])
        verbose = sys.argv[4]

    mainFunction(instanceFilename, assignmentFilename, verbose, timeLimit)
