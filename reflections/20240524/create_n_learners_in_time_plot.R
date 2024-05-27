#!/bin/env Rscript
n_learners_in_time <- readr::read_csv(
  "n_learners_in_time.csv"
  )

ggplot2::ggplot(
  data = n_learners_in_time,
  mapping = ggplot2::aes(x = time, y = n_learners)
) + ggplot2::geom_line(
) + ggplot2::geom_rect(mapping = ggplot2::aes(
    xmin = readr::parse_time("9:00"), 
    xmax = readr::parse_time("10:20"), 
    ymin = 0, 
    ymax = 8
  ), fill = ggplot2::alpha("orange", 0.01)
) + ggplot2::geom_rect(mapping = ggplot2::aes(
  xmin = readr::parse_time("15:30"), 
  xmax = readr::parse_time("16:00"), 
  ymin = 0, 
  ymax = 8
), fill = ggplot2::alpha("orange", 0.01)
) + ggplot2::geom_rect(mapping = ggplot2::aes(
  xmin = readr::parse_time("12:06"), 
  xmax = readr::parse_time("13:00"), 
  ymin = 0, 
  ymax = 8
), fill = ggplot2::alpha("blue", 0.01)
) + ggplot2::theme(text = ggplot2::element_text(size = 20)) +
  ggplot2::labs(
  title = "Number of learners in time",
  caption = "Intermediate Bianca, 2024-05-25, orange = me, blue = break"
)
ggplot2::ggsave("n_learners_in_time.png", width = 7, height = 7)
