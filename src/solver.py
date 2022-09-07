#! usr/bin/python
# -*- coding: utf-8 -*-

from mip import *
import problem as pb


def solve(data: pb.Data, maxTime: int, verbose: bool) -> pb.Solution:

    model = Model(sense=MINIMIZE, solver_name=GRB)
    model.verbose = int(verbose)

    # Variables:
    x = [
        [model.add_var(var_type=BINARY) for i in range(data.nbMachines)]
        for j in range(data.nbProcess)
    ]

    y = [
        [model.add_var(var_type=BINARY) for i in range(data.nbLocations)]
        for j in range(data.nbServices)
    ]

    z = [
        [model.add_var(var_type=BINARY) for i in range(data.nbNeighborhoods)]
        for j in range(data.nbServices)
    ]

    lc_1 = [
        [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbResources)]
        for j in range(data.nbMachines)
    ]

    lc_2 = [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbResources)]

    a = [
        [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbResources)]
        for j in range(data.nbMachines)
    ]

    bc_1 = [
        [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbMachines)]
        for j in range(data.nbBalanceTriples)
    ]

    bc_2 = [
        model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbBalanceTriples)
    ]

    smc_1 = [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbServices)]

    smc_2 = model.add_var(var_type=CONTINUOUS, lb=0)

    mmc_1 = [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbProcess)]

    # Constraints:
    # Assignment
    for p in range(data.nbProcess):
        model += xsum(x[p][m] for m in range(data.nbMachines)) == 1, "Assignment"

    # Capacity
    for r in range(data.nbResources):
        for m in range(data.nbMachines):
            model += (
                xsum(data.processReq[p][r] * x[p][m] for p in range(data.nbProcess))
                <= data.hardResCapacities[m][r]
            ), "Capacity"

    # Conflict
    for s in range(data.nbServices):
        for m in range(data.nbMachines):
            model += (
                xsum(
                    x[p][m]
                    for p in range(data.nbServices)
                    if data.servicesProcess[p] == s
                )
                <= 1
            ), "Conflict"

    # Spread 1
    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            model += (
                xsum(
                    xsum(
                        x[p][m]
                        for p in range(data.nbServices)
                        if data.servicesProcess[p] == s
                    )
                    for m in range(data.nbLocations)
                    if data.locations[m] == l
                )
                <= data.nbLocations * data.nbServices * y[s][l]
            ), "Spread 1"

    # Spread 2
    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            model += (
                xsum(
                    xsum(
                        x[p][m]
                        for p in range(data.nbServices)
                        if data.servicesProcess[p] == s
                    )
                    for m in range(data.nbLocations)
                    if data.locations[m] == l
                )
                >= y[s][l],
                "Spread 2",
            )

    # Spread 3
    for s in range(data.nbServices):
        model += (
            xsum(y[s][l] for l in range(data.nbLocations)) >= data.spreadMin[s],
            "Spread 3",
        )

    # Transient
    for r in range(data.nbResources):
        if data.transientStatus[r] == 1:
            for m in range(data.nbMachines):
                model += (
                    xsum(
                        data.processReq[p][r]
                        for p in range(data.nbProcess)
                        if data.initialAssignment[p] == m
                    )
                    + xsum(
                        data.processReq[p][r] * x[p][m]
                        for p in range(data.nbProcess)
                        if data.initialAssignment[p] != m
                    )
                    <= data.hardResCapacities[m][r]
                ), "Transient"

    # Dependency 1
    for s1 in range(data.nbServices):
        for s2 in data.dependencies[s1]:
            for n in range(data.nbNeighborhoods):
                model += z[s1][n] <= z[s2][n], "Dependency 1"

    # Dependency 2
    for s in range(data.nbServices):
        for n in range(data.nbNeighborhoods):
            for p in range(data.nbProcess):
                if data.servicesProcess[p] == s:
                    for m in range(data.nbMachines):
                        if data.neighborhoods[m] == n:
                            model += (
                                x[p][m] <= z[s][n],
                                "Dependency 2",
                            )

    # Dependency 3
    for s in range(data.nbServices):
        for n in range(data.nbNeighborhoods):
            model += (
                xsum(
                    xsum(
                        x[p][m]
                        for m in range(data.nbMachines)
                        if data.neighborhoods[m] == n
                    )
                    for p in range(data.nbProcess)
                    if data.servicesProcess[p] == s
                )
                >= z[s][n]
            ), "Dependency 3"

    # Objective constraints:
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
            (1 - x[p][data.initialAssignment[p]])
            for p in range(data.nbServices)
            if data.servicesProcess[p] == s
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

    # Objective
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
    model.write("../model.lp")
    model.read("../model.lp")
    if verbose:
        print(
            "model has {} vars, {} constraints and {} nzs".format(
                model.num_cols, model.num_rows, model.num_nz
            )
        )

    # Solve
    status = model.optimize(max_seconds=maxTime)
    if status == OptimizationStatus.OPTIMAL:
        print("optimal solution cost {} found".format(model.objective_value))
    elif status == OptimizationStatus.FEASIBLE:
        print(
            "sol.cost {} found, best possible: {}".format(
                model.objective_value, model.objective_bound
            )
        )
    elif status == OptimizationStatus.NO_SOLUTION_FOUND:
        print(
            "no feasible solution found, lower bound is: {}".format(
                model.objective_bound
            )
        )

    assignment = []
    for p in range(data.nbProcess):
        for m in range(data.nbMachines):
            if x[p][m].x > 0.9:
                assignment.append(m)

    return pb.Solution(assignment, model.objective_value)
