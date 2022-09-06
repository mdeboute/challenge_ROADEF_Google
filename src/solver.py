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
    for p in range(data.nbProcess):
        model += xsum(x[p][m] for m in range(data.nbMachines)) == 1

    for r in range(data.nbResources):
        for m in range(data.nbMachines):
            model += (
                xsum(data.processReq[p][r] * x[p][m] for m in range(data.nbMachines))
                <= data.hardResCapacities[m][r]
            )

    return
