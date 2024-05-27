#!/bin/env Rscript
n_learners_in_time <- readr::read_csv(
  "20240524_n_learners_in_time.csv"
  )

# events <- tibble::tibble(
#   xmin = readr::parse_time(c("9:00", "10:00", "10:15", "11:00", "12:06", "14:00", "15:00", "15:30")),
#   xmax = readr::parse_time(c("10:00","10:15", "10:35", "11:15", "13:00", "14:15", "15:15", "16:00")),
#   event =                 c("orange", "blue", "orange", "blue",  "blue", "blue",  "blue",   "orange")
# )

events <- tibble::tribble(
  ~xmin, ~xmax, ~event,
  "9:00", "10:00", "orange",
  "10:00", "10:15", "blue",
  "10:15", "10:35", "orange",
  "11:00", "11:15", "blue",
  "12:06", "13:00", "blue",
  "14:00", "14:15", "blue",
  "15:00", "15:15", "blue",
  "15:30", "16:00", "orange"
)
events$xmin <- readr::parse_time(events$xmin)
events$xmax <- readr::parse_time(events$xmax)

ggplot2::ggplot(
  data = n_learners_in_time,
  mapping = ggplot2::aes(x = time, y = n_learners)
) + ggplot2::geom_line() + 
  ggplot2::geom_point(size = 5) + 
  ggplot2::scale_y_continuous(name = "Number of learners") +
  ggplot2::scale_x_time(
    name = "Time", 
    breaks = seq(
      from = readr::parse_time("9:00"), 
      to = readr::parse_time("16:00"), 
      length.out = 8
    )
  ) +
  ggplot2::theme(text = ggplot2::element_text(size = 20)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::labs(
    title = "Number of learners in time",
    caption = "Intermediate Bianca, 2024-05-25, orange = R, blue = break"
  ) + ggplot2::geom_rect(
    data = events,
    inherit.aes = FALSE, 
    ggplot2::aes(
      xmin = xmin,
      ymin = 0,
      xmax = xmax,
      ymax = 8
    ),
    fill = ggplot2::alpha(events$event, 0.5)
  )

ggplot2::ggsave("20240524_n_learners_in_time.png", width = 7, height = 7)
