# Replicate jobs

!!!- info "Learning objectives"

    - Practice using the UPPMAX documentation
    - I can schedule a job that has replicates

???- question "For teachers"

    Teaching goals are:

    - Learners have practiced using the UPPMAX documentation
    - Learners have scheduled a job that has replicates

    Lesson plan:

    ```mermaid
    gantt
      title Replicate jobs
      dateFormat X
      axisFormat %s
      section First hour
      Course introduction: done, course_intro, 0, 10s
      Prior : intro, after course_intro, 5s
      Present: theory_1, after intro, 5s
      Challenge: crit, exercise_1, after theory_1, 40s
      Break: crit, milestone, after exercise_1
      section Second hour
      Challenge: crit, exercise_2, 0, 10s
      Feedback: feedback_2, after exercise_2, 10s
      SLURM: done, slurm, after feedback_2, 25s
      Break: done, milestone, after slurm
    ```

    Prior questions:

    - You have a simulation that uses randomness. You want to run
      it 1000 times. How would you do that?


## Why?

You have a simulation that uses randomness. 
You want to run 1000 replicates of it. How would you do that?

## `run_simulation.sh`

```bash
#!/bin/bash
echo $((1 + ($RANDOM % 6))) > result_$1.txt
```

## `submit_runs.sh`

```bash
#!/bin/bash
#SBATCH --array=0-10
./run_simulation.sh ${SLURM_ARRAY_TASK_ID}
```

## Exercises

