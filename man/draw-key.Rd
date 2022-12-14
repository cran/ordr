% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biplot-key.r
\name{draw-key}
\alias{draw-key}
\alias{draw_key_line}
\alias{draw_key_crosslines}
\alias{draw_key_crosspoint}
\title{Biplot key drawing functions}
\usage{
draw_key_line(data, params, size)

draw_key_crosslines(data, params, size)

draw_key_crosspoint(data, params, size)
}
\arguments{
\item{data}{A single row data frame containing the scaled aesthetics to
display in this key}

\item{params}{A list of additional parameters supplied to the geom.}

\item{size}{Width and height of key in mm.}
}
\value{
A grid grob.
}
\description{
These key drawing functions supplement those built into
\strong{\link[ggplot2:ggplot2-package]{ggplot2}} for producing legends suitable to biplots.
}
\details{
\code{draw_key_line()} is a horizontal counterpart to \code{\link[ggplot2:draw_key]{ggplot2::draw_key_vline()}}.
\code{draw_key_crosslines()} superimposes these two keys, and
\code{draw_key_crosspoint()} additionally superimposes an oversized
\code{\link[ggplot2:draw_key]{ggplot2::draw_key_point()}}.
}
\examples{
# scaled PCA of Anderson iris data with ranges and confidence intervals
iris[, -5] \%>\%
  prcomp(scale = TRUE) \%>\%
  as_tbl_ord() \%>\%
  confer_inertia(1) \%>\%
  augment_ord() \%>\%
  mutate_rows(species = iris$Species) \%>\%
  ggbiplot(aes(color = species)) +
  theme_bw() +
  scale_color_brewer(type = "qual", palette = 2) +
  geom_rows_lineranges(fun.data = mean_sdl, linewidth = .75) +
  geom_rows_point(alpha = .5) +
  geom_cols_vector(color = "#444444") +
  geom_cols_text_radiate(aes(label = name), color = "#444444", size = 3) +
  ggtitle(
    "Row-principal PCA biplot of Anderson iris data",
    "Ranges 2 sample standard deviations from centroids"
  )
}
\seealso{
\link[ggplot2:draw_key]{ggplot2::draw_key} for key glyphs installed with \strong{ggplot2}.
}
