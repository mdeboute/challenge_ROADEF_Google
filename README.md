# challenge_ROADEF_Google

This project is the ROADEF/EURO challenge of 2012. It is a competition jointly organized by the French Operational Research and Decision Aid society ([ROADEF](https://www.roadef.org/societe-francaise-recherche-operationnelle-aide-decision)) and the European Operational Research society (EURO). This competition aims at developing a solution method for an industrial optimization problem proposed by an industrial partner. In 2012, the problem was proposed by Google. This problem was a large-scale machine reassignment problem. This project was supervised by Mr. Eduardo Uchoa.

[Here](./Problem-definition-Google-ROADEF-challenge.pdf) is the wording of the problem.

Solve one instance:

```shell
Usage: python3 src/run.py <model_file> <assignment_file> <time_limit (s)> <verbose (int)>
```

Solve all instances of a directory:

```shell
Usage: ./benchmark.sh <data_dir> <time_limit (s)> <verbose (int)>
```

## Beamer

[Here](./beamer.pdf) is a small presentation of our model.
