#! usr/bin/python
# -*- coding: utf-8 -*-

from mip import *
import problem as pb
import numpy as np
import time


def solve(data: pb.Data, maxTime: int, verbose: bool) -> pb.Solution:

    model = Model(name="GoogleChallenge", sense=MINIMIZE, solver_name=GRB)
    model.verbose = int(verbose)

    ### Decision Variables ###

    # main variable indicating if a process p is assigned to the machine m
    x = np.array(
        [
            [
                model.add_var(var_type=BINARY, name="x(" + str(p) + "," + str(m) + ")")
                for m in range(data.nbMachines)
            ]
            for p in range(data.nbProcess)
        ]
    )

    # auxiliary variable indicating if service s has at least one process in location l
    y = np.array(
        [
            [
                model.add_var(var_type=BINARY, name="y(" + str(s) + "," + str(l) + ")")
                for l in range(data.nbLocations)
            ]
            for s in range(data.nbServices)
        ]
    )

    # auxiliary variable indicating if service s has at least one process in neighborhood n
    z = np.array(
        [
            [
                model.add_var(var_type=BINARY, name="z(" + str(s) + "," + str(n) + ")")
                for n in range(data.nbNeighborhoods)
            ]
            for s in range(data.nbServices)
        ]
    )

    ### Cost Variables ###

    # auxiliary variable to compute the used capacity beyond the safety capacity 
    # of resource r on machine m
    lc_1 = np.array(
        [
            [
                model.add_var(
                    var_type=CONTINUOUS,
                    lb=0,
                    name="lc_1(" + str(m) + "," + str(r) + ")",
                )
                for r in range(data.nbResources)
            ]
            for m in range(data.nbMachines)
        ]
    )

    # auxiliary variable to compute the load cost of resource r on all machines
    lc_2 = np.array(
        [
            model.add_var(var_type=CONTINUOUS, lb=0, name="lc_2(" + str(r) + ")")
            for r in range(data.nbResources)
        ]
    )

    # auxiliary variable to compute the free capacity of resource r on machine m
    a = np.array(
        [
            [
                model.add_var(
                    var_type=CONTINUOUS, lb=0, name="a(" + str(m) + "," + str(r) + ")"
                )
                for r in range(data.nbResources)
            ]
            for m in range(data.nbMachines)
        ]
    )

    # auxiliary variable to compute the balance of resources on machine m for triplet b
    bc_1 = np.array(
        [
            [
                model.add_var(
                    var_type=CONTINUOUS,
                    lb=0,
                    name="bc_1(" + str(b) + "," + str(m) + ")",
                )
                for m in range(data.nbMachines)
            ]
            for b in range(data.nbBalanceTriples)
        ]
    )

    # auxiliary variable to compute the balance cost of triplet b
    bc_2 = np.array(
        [
            model.add_var(var_type=CONTINUOUS, lb=0, name="bc_2(" + str(b) + ")")
            for b in range(data.nbBalanceTriples)
        ]
    )

    # auxiliary variable to compute the amount of moved processes of service s
    smc_1 = np.array(
        [
            model.add_var(var_type=CONTINUOUS, lb=0, name="smc_1(" + str(s) + ")")
            for s in range(data.nbServices)
        ]
    )

    # auxiliary variable to compute the service move cost
    smc_2 = model.add_var(var_type=CONTINUOUS, lb=0, name="smc_2")

    # auxiliary variable to compute the machine move cost of process p
    mmc_1 = np.array(
        [
            model.add_var(var_type=CONTINUOUS, lb=0, name="mmc_1(" + str(p) + ")")
            for p in range(data.nbProcess)
        ]
    )

    # auxiliary variable to compute the machine move cost
    mmc_2 = model.add_var(var_type=CONTINUOUS, lb=0, name="mmc_2")


    ### Constraints ###
    # Sets definition:
    initialProcessesInMachine = np.frompyfunc(list, 0, 1)(
        np.empty(data.nbMachines, dtype=object)
    )
    initialProcessesNotInMachine = np.frompyfunc(list, 0, 1)(
        np.empty(data.nbMachines, dtype=object)
    )
    processesOfService = np.frompyfunc(list, 0, 1)(
        np.empty(data.nbServices, dtype=object)
    )
    machinesOfLocation = np.frompyfunc(list, 0, 1)(
        np.empty(data.nbLocations, dtype=object)
    )
    machinesOfNeighbourhood = np.frompyfunc(list, 0, 1)(
        np.empty(data.nbNeighborhoods, dtype=object)
    )

    for m in range(data.nbMachines):
        for p in range(data.nbProcess):
            if data.initialAssignment[p] == m:
                initialProcessesInMachine[m].append(p)
            else:
                initialProcessesNotInMachine[m].append(p)

        for l in range(data.nbLocations):
            if data.locations[m] == l:
                machinesOfLocation[l].append(m)

        for n in range(data.nbNeighborhoods):
            if data.neighborhoods[m] == n:
                machinesOfNeighbourhood[n].append(m)

    for p in range(data.nbProcess):
        for s in range(data.nbServices):
            if data.servicesProcess[p] == s:
                processesOfService[s].append(p)

    # Assignment
    for p in range(data.nbProcess):
        model += xsum(
            x[p][m] for m in range(data.nbMachines)
        ) == 1, "Assignment_" + str(p)

    # Capacity
    for r in range(data.nbResources):
        for m in range(data.nbMachines):
            model += (
                xsum(data.processReq[p][r] * x[p][m] for p in range(data.nbProcess))
                <= data.hardResCapacities[m][r]
            ), "Capacity_" + str(r) + "_" + str(m)

    # Conflict
    for s in range(data.nbServices):
        for m in range(data.nbMachines):
            model += (
                xsum(x[p][m] for p in processesOfService[s]) <= 1
            ), "Conflict_" + str(s) + "_" + str(m)

    # Spread 1
    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            for p in processesOfService[s]:
                for m in machinesOfLocation[l]:
                    model += (x[p][m] <= y[s][l], "Spread_1_" + str(l) + "_" + str(s))

    # Spread 2
    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            model += (
                xsum(
                    xsum(x[p][m] for p in processesOfService[s])
                    for m in machinesOfLocation[l]
                )
                >= y[s][l],
                "Spread_2_" + str(l) + "_" + str(s),
            )

    # Spread 3
    for s in range(data.nbServices):
        model += (
            xsum(y[s][l] for l in range(data.nbLocations)) >= data.spreadMin[s],
            "Spread_3_" + str(s),
        )

    # Transient
    for r in range(data.nbResources):
        if data.transientStatus[r] == 1:
            for m in range(data.nbMachines):
                model += (
                    xsum(data.processReq[p][r] for p in initialProcessesInMachine[m])
                    + xsum(
                        data.processReq[p][r] * x[p][m]
                        for p in initialProcessesNotInMachine[m]
                    )
                    <= data.hardResCapacities[m][r]
                ), "Transient" + str(r) + "_" + str(m)

    # Dependency 1
    for s1 in range(data.nbServices):
        for s2 in data.dependencies[s1]:
            for n in range(data.nbNeighborhoods):
                model += z[s1][n] <= z[s2][n], "Dependency_1_" + str(s1) + "_" + str(
                    s2
                ) + "_" + str(n)

    # Dependency 2
    for s in range(data.nbServices):
        for n in range(data.nbNeighborhoods):
            for p in processesOfService[s]:
                for m in machinesOfNeighbourhood[n]:
                    model += (
                        (x[p][m] <= z[s][n]),
                        "Dependency_2_"
                        + str(s)
                        + "_"
                        + str(n)
                        + "_"
                        + str(p)
                        + "_"
                        + str(m),
                    )

    # Dependency 3
    for s in range(data.nbServices):
        for n in range(data.nbNeighborhoods):
            model += (
                xsum(
                    xsum(x[p][m] for m in machinesOfNeighbourhood[n])
                    for p in processesOfService[s]
                )
                >= z[s][n]
            ), "Dependency_3_" + str(s) + "_" + str(n)

    ### Cost computation constraints ###
    # Load cost
    for r in range(data.nbResources):
        for m in range(data.nbMachines):
            model += (
                lc_1[m][r]
                >= xsum(data.processReq[p][r] * x[p][m] for p in range(data.nbProcess))
                - data.softResCapacities[m][r]
            )

    for r in range(data.nbResources):
        lc_2[r] = xsum(lc_1[m][r] for m in range(data.nbMachines))

    # Balance cost
    for m in range(data.nbMachines):
        for r in range(data.nbResources):
            a[m][r] = data.hardResCapacities[m][r] - xsum(
                data.processReq[p][r] * x[p][m] for p in range(data.nbProcess)
            )

    for b in range(data.nbBalanceTriples):
        for m in range(data.nbMachines):
            model += (
                bc_1[b][m]
                >= data.balanceTriples[b].target
                * a[m][data.balanceTriples[b].resource1]
                - a[m][data.balanceTriples[b].resource2]
            )

    for b in range(data.nbBalanceTriples):
        bc_2[b] = xsum(bc_1[b][m] for m in range(data.nbMachines))

    # Process move cost
    pmc = xsum(
        (1 - x[p][data.initialAssignment[p]]) * data.processMoveCost[p]
        for p in range(data.nbProcess)
    )

    # Service move cost
    for s in range(data.nbServices):
        smc_1[s] = xsum(
            (1 - x[p][data.initialAssignment[p]]) for p in processesOfService[s]
        )

    for s in range(data.nbServices):
        model += smc_2 >= smc_1[s]

    # Machine move cost
    for p in range(data.nbProcess):
        mmc_1[p] = xsum(
            x[p][m] * data.machineMoveCosts[data.initialAssignment[p]][m]
            for m in range(data.nbMachines)
        )

    mmc_2 = xsum(mmc_1[p] for p in range(data.nbProcess))

    ### Objective ###
    model.objective = minimize(
        xsum(data.weightLoadCost[r] * lc_2[r] for r in range(data.nbResources))
        + xsum(
            data.balanceTriples[b].weight * bc_2[b]
            for b in range(data.nbBalanceTriples)
        )
        + data.processMoveWeight * pmc
        + data.serviceMoveWeight * smc_2
        + data.machineMoveWeight * mmc_2
    )

    # Save model
    model.write("model.lp")

    # Limitation of the number of processors
    model.threads = 1
    model.max_seconds = maxTime
    model.max_mip_gap = 1e-8
    model.max_mip_gap_abs = 1
    model.solver.set_verbose(verbose)

    # Launching the stopwatch
    start = time.perf_counter()

    # Model resolution
    status = model.optimize()

    # Stopping the stopwatch and resolution time calculation
    runtime = time.perf_counter() - start

    print("\n----------------------------------")
    if status == OptimizationStatus.OPTIMAL:
        print("Optimization status: OPTIMAL")
    elif status == OptimizationStatus.FEASIBLE:
        print(
            "Optimization status: TIME LIMIT and A FEASIBLE SOLUTION COMPUTED (WITHOUT OPTIMALITY GUARANTEE)"
        )
    elif status == OptimizationStatus.NO_SOLUTION_FOUND:
        print("Optimization status: TIME LIMIT and NO SOLUTION COMPUTED")
    elif (
        status == OptimizationStatus.INFEASIBLE
        or status == OptimizationStatus.INT_INFEASIBLE
    ):
        print("Optimization status: INFEASIBLE")
    elif status == OptimizationStatus.UNBOUNDED:
        print("Optimization status: UNBOUNDED")

    print("Solution time (s) : ", runtime)
    print("----------------------------------")

    assignment = np.empty(data.nbProcess, dtype=int)
    val = -1
    # If the model has been solved to optimality or if a solution has been found within the given time limit
    if model.num_solutions > 0:
        print(
            "Objective function value of the computed solution: ",
            round(model.objective_value),
        )
        val = round(model.objective_value)
        print(
            "Best lower bound on the optimal objective function value: ",
            model.objective_bound,
        )
        for p in range(data.nbProcess):
            for m in range(data.nbMachines):
                # 0.9 is used since the value of binary variables can sometimes be 0.99999 or 0.00001
                if x[p][m].x > 0.9 and x[p][m].name.startswith("x"):
                    if verbose:
                        if data.initialAssignment[p] != m:
                            print("\t Process ", p, " moved to machine ", m)
                        else:
                            print("\t Process ", p, " remained in machine ", m)
                    assignment[p] = m
    else:
        print("No solution computed")
    print("----------------------------------\n")

    return pb.Solution(assignment, val)
