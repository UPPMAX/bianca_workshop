#!/bin/bash
job_id_a=$(sbatch -A sens2023598 do_a.sh | cut -d " " -f 4)
job_id_b=$(sbatch -A sens2023598 do_b.sh | cut -d " " -f 4)
sbatch -A sens2023598 do_c.sh --dependency=afterok:${job_id_a}:${job_id_b}
