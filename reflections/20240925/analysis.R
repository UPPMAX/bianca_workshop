# I can transfer files to/from Rackham using FileZilla,0,1
# I can transfer files to/from Rackham using FileZilla,2,1
# I can transfer files to/from Rackham using FileZilla,3,3
# I can transfer files to/from Rackham using FileZilla,4,3
# I can transfer files to/from Rackham using FileZilla,5,9
confidences_transfer_r <- c(rep(0, 1), rep(2, 1), rep(3, 3), rep(4,3), rep(5, 9))
testthat::expect_equal(length(confidences_transfer_r), 17)
testthat::expect_equal(mean(confidences_transfer_r), 4.0)

# I can start an interactive session,4,2
# I can start an interactive session,5,15
# I can start an interactive session,4,2
# I can start an interactive session,5,15
# Average = ((4 * 2) + (5 * 15)) / 17 = 4.88
confidences_interactive_r <- c(rep(4, 2), rep(5, 15))
testthat::expect_equal(length(confidences_interactive_r), 17)
testthat::expect_equal(mean(confidences_interactive_r), 4.88, tolerance = 0.01)

# I can schedule a job,2,1
# I can schedule a job,3,1
# I can schedule a job,4,3
# I can schedule a job,5,12
# Average = ((2 * 1) + (3 * 1) + (4 * 3) + (5 * 12) ) / 17 = 4.53
confidences_schedule_r <- c(rep(2, 1), rep(3, 1), rep(4, 3), rep(5, 12))
testthat::expect_equal(length(confidences_schedule_r), 17)
testthat::expect_equal(mean(confidences_schedule_r), 4.53, tolerance = 0.01)

## 2024-09-25: Bianca Intro

# I can transfer files to/from Bianca using FileZilla,1,1
# I can transfer files to/from Bianca using FileZilla,2,2
# I can transfer files to/from Bianca using FileZilla,3,1
# I can transfer files to/from Bianca using FileZilla,4,1
# Average = ((1 * 1) + (2 * 2) + (3 * 1) + (4 * 1)) / 5 = 2.4
confidences_transfer_nr <- c(rep(1, 1), rep(2, 2), rep(3, 1), rep(4,1))
testthat::expect_equal(length(confidences_transfer_nr), 5)
testthat::expect_equal(mean(confidences_transfer_nr), 2.4)

# I can start an interactive session,4,3
# I can start an interactive session,5,2
# Average = ((4 * 3) + (5 * 2)) / 5 = 4.4
confidences_interactive_nr <- c(rep(4, 3), rep(5, 2))
testthat::expect_equal(length(confidences_interactive_nr), 5)
testthat::expect_equal(mean(confidences_interactive_nr), 4.4, tolerance = 0.01)

# I can submit jobs to the scheduler,3,1
# I can submit jobs to the scheduler,4,2
# I can submit jobs to the scheduler,5,2
# Average = ((3 * 1) + (4 * 2) + (5 * 2)) / 5 = 4.2
confidences_schedule_nr <- c(rep(3, 1), rep(4, 2), rep(5, 2))
testthat::expect_equal(length(confidences_schedule_nr), 5)
testthat::expect_equal(mean(confidences_schedule_nr), 4.2, tolerance = 0.01)

ks_test <- ks.test(rpois(n=20, lambda=5), rpois(n=20, lambda=5))
testthat::expect_true(ks_test$p.value >= 0.05) # (usually) same distribution

ks_transfer <- ks.test(confidences_transfer_r, confidences_transfer_nr)
testthat::expect_true(ks_transfer$p.value >= 0.05) # Same distribution

ks_interactive <- ks.test(confidences_interactive_r, confidences_interactive_nr)
testthat::expect_true(ks_interactive$p.value >= 0.05) # Same distribution

ks_schedule <- ks.test(confidences_schedule_r, confidences_schedule_nr)
testthat::expect_true(ks_schedule$p.value >= 0.05) # Same distribution
