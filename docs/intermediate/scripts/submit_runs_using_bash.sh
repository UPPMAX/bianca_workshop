#!/bin/bash
for i in {0..10}
do
  sbatch run_simulation.sh ${i}
done
