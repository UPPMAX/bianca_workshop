#!/bin/bash
job_id_a=$(sbatch -A staff do_a.sh | cut -d " " -f 4)
job_id_b=$(sbatch -A staff do_b.sh | cut -d " " -f 4)
sbatch -A staff --dependency=afterok:${job_id_a},${job_id_b} do_c.sh
