#! usr/bin/python
# -*- coding: utf-8 -*-

from mip import *
import problem as pb


def solve(data: pb.Data, maxTime: int, verbose: bool) -> pb.Solution:

    model = Model(sense=MINIMIZE, solver_name=GRB)

    # Variables
    x = [
        [model.add_var(var_type=BINARY) for i in range(data.nbProcess)]
        for j in range(data.nbMachines)
    ]

    y = [
        [model.add_var(var_type=BINARY) for i in range(data.nbServices)]
        for j in range(data.nbLocations)
    ]

    z = [
        [model.add_var(var_type=BINARY) for i in range(data.nbServices)]
        for j in range(data.nbNeighborhoods)
    ]

    lc = [
        [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbServices)]
        for j in range(data.nbLocations)
    ]

    bc = [
        [model.add_var(var_type=CONTINUOUS, lb=0) for i in range(data.nbBalanceTriples)]
        for j in range(data.nbMachines)
    ]

    smc = model.add_var(var_type=CONTINUOUS, lb=0)

    # Constraints

    # Assignment
    for p in range(data.nbProcess):
        model += xsum(x[p][m] for m in range(data.nbMachines)) == 1

    # Capacity
    for r in range(data.nbResources):
        for m in range(data.nbMachines):
            model += (
                xsum(data.processReq[p][r] * x[p][m] for p in range(data.nbProcess))
                <= data.hardResCapacities[m][r]
            )

    # Conflict
    for s in range(data.nbServices):
        for m in range(data.nbMachines):
            model += xsum(x[p][m] for p in s) <= 1

    # Spread
    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            model += (
                xsum(xsum(x[p][m] for p in s) for m in l)
                <= data.nbLocations * data.nbServices * y[s][l]
            )

    for l in range(data.nbLocations):
        for s in range(data.nbServices):
            model += xsum(xsum(x[p][m] for p in s) for m in l) >= y[s][l]

    for s in range(data.nbServices):
        model += xsum(y[s][l] for l in range(data.nbLocations)) >= data.spreadMin[s]

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
                )

    # Dependency 1

    return
