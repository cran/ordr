## ----knitr options, echo=FALSE, message=FALSE, warning=FALSE------------------
knitr::opts_chunk$set(
  message = FALSE,
  fig.height = 6, fig.align = "center"
)

## ----glass pca----------------------------------------------------------------
data(glass, package = "ordr")
x <- scale(glass[, c("SiO2", "Al2O3", "FeO", "MgO", "CaO")],
           center = TRUE, scale = TRUE)
s <- svd(x)
r <- length(s$d)

## ----inertia------------------------------------------------------------------
# inertia of the (scaled) data
sum(x^2)
# inertia of the case and variable factors
sum(s$u^2)
sum(s$v^2)
# inertia of the diagonal factor
sum(s$d^2)

## ----case geometry, fig.width=5.5---------------------------------------------
# distances between cases
x.dist <- dist(x)
# distances between cases (principal coordinates)
s.dist <- dist(s$u[, 1:2] %*% diag(s$d[1:2]))
# scatterplot
plot(
  x = as.vector(x.dist),
  y = as.vector(s.dist),
  xlim = c(0, 10), ylim = c(0, 10),
  asp = 1, pch = 19, cex = .5,
  xlab = "Case distances in centered and scaled data",
  ylab = "Case point distances in planar biplot"
)
lines(x = c(0, 10), y = c(0, 10))

## ----variable geometry, fig.width=5.5-----------------------------------------
# correlations between variables
x.cor <- cor(x)
# magnitudes of variable vectors
s.len <- apply(s$v[, 1:2] %*% diag(s$d[1:2]), 1, norm, "2")
# cosines between variables (principal coordinates)
s.cor <- (s$v[, 1:2] / s.len) %*% diag(s$d[1:2]^2) %*% t(s$v[, 1:2] / s.len)
# scatterplot
plot(
  x = as.vector(x.cor[lower.tri(x.cor)]),
  y = as.vector(s.cor[lower.tri(s.cor)]),
  xlim = c(-1, 1), ylim = c(-1, 1),
  asp = 1, pch = 19, cex = .5,
  xlab = "Variable correlations in centered and scaled data",
  ylab = "Variable vector cosines in planar biplot"
)
lines(x = c(-1, 1), y = c(-1, 1))

## ----multidimensional scaling of cities, fig.width=5.5------------------------
d <- as.matrix(UScitiesD)
cent <- diag(1, nrow(d)) - matrix(1/nrow(d), nrow(d), nrow(d))
d.cent <- -.5 * cent %*% (d^2) %*% cent
d.cmds <- svd(d.cent)
d.coord <- d.cmds$u[, 1:2] %*% diag(sqrt(d.cmds$d[1:2]))
# scatterplot
plot(
  x = as.vector(UScitiesD),
  y = as.vector(dist(d.coord)),
  xlim = c(0, max(UScitiesD)), ylim = c(0, max(UScitiesD)),
  asp = 1, pch = 19, cex = .5,
  xlab = "City road distances",
  ylab = "Point distances in planar CMDS"
)
lines(x = c(0, max(UScitiesD)), y = c(0, max(UScitiesD)))

## ----multidimensional scaling of cities plot, fig.height=5, fig.width=7-------
plot(
  d.coord, pch = NA, asp = 1,
  xlab = "First principal coordinate", ylab = "Second principal coordinate"
)
text(d.coord, labels = rownames(d), cex = .9)

## ----multidimensional scaling of glass, fig.width=5.5-------------------------
# covariances and standard deviations
c <- cov(x)
s <- diag(sqrt(diag(c)))
# eigendecomposition of covariance matrix
c.eigen <- eigen(c)
# artificial coordinates
c.coord <- c.eigen$vectors[, 1:2] %*% diag(sqrt(c.eigen$values[1:2]))
# scatterplot
c.inner <- c.coord %*% t(c.coord)
plot(
  x = as.vector(c[lower.tri(c)]),
  y = as.vector(c.inner[lower.tri(c.inner)]),
  xlim = range(c[lower.tri(c)]), ylim = range(c[lower.tri(c)]),
  asp = 1, pch = 19, cex = .5,
  xlab = "Measurement covariances in unscaled data",
  ylab = "Vector inner products in planar CMDS"
)
lines(x = range(c[lower.tri(c)]), y = range(c[lower.tri(c)]))

## ----multidimensional scaling of glass plot, fig.width=7----------------------
c <- cor(glass[, c("SiO2", "Al2O3", "FeO", "MgO", "CaO")])
c.eigen <- eigen(c)
c.coord <- c.eigen$vectors[, 1:2] %*% diag(sqrt(c.eigen$values[1:2]))
plot(
  c.coord, pch = NA, asp = 1,
  xlab = "First principal coordinate", ylab = "Second principal coordinate"
)
arrows(0, 0, c.coord[, 1L], c.coord[, 2L])
text(c.coord, labels = rownames(c), cex = .9)

## ----big guns-----------------------------------------------------------------
library(ordr)
library(dplyr)

## ----QS top university rankings-----------------------------------------------
data(qswur_usa, package = "ordr")
print(qswur_usa)

## ----calibrate rankings-------------------------------------------------------
qswur_usa %>%
  filter(year == 2020L) %>%
  select(institution, starts_with("rk_")) %>%
  mutate_at(
    vars(starts_with("rk_")),
    ~ match(., sort(unique(as.numeric(.))))
  ) %>%
  filter_at(vars(starts_with("rk_")), ~ ! is.na(.)) ->
  qswur_usa2020
print(qswur_usa2020)

## ----Kendall rank correlations, fig.width=7-----------------------------------
corr <- cor(select(qswur_usa2020, starts_with("rk_")), method = "kendall")
heatmap(corr, scale = "none")

## ----multidimensional scaling of variables, fig.width=7, eval=FALSE, echo=FALSE----
# corr.eigen <- eigen(corr)
# corr.coord <- corr.eigen$vectors %*% diag(sqrt(corr.eigen$values))
# plot(corr.coord, pch = NA, asp = 1, xlab = "", ylab = "")
# arrows(0, 0, corr.coord[, 1], corr.coord[, 2])
# text(corr.coord, labels = rownames(corr))

## ----fig.width=7--------------------------------------------------------------
eigen_ord(corr) %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  mutate_rows(metric = stringr::str_remove(name, "rk_")) %>%
  confer_inertia(1) ->
  c_eigen
c_eigen %>%
  ggbiplot() +
  theme_minimal() +
  geom_unit_circle() +
  geom_rows_vector(aes(label = metric)) +
  scale_x_reverse(expand = expansion(add = .4)) +
  scale_y_continuous(expand = expansion(add = .3)) +
  ggtitle("Kendall correlations between university rankings",
          "CMDS correlation monoplot")

## -----------------------------------------------------------------------------
c_eigen %>%
  fortify(.matrix = "rows") %>%
  select(-name, -.matrix)

## ----fig.width=7--------------------------------------------------------------
c_eigen %>%
  ggbiplot(aes(x = 2, y = 3)) +
  theme_minimal() +
  geom_unit_circle() +
  geom_rows_vector(aes(label = metric)) +
  scale_x_continuous(expand = expansion(add = .5)) +
  scale_y_continuous(expand = expansion(add = .5)) +
  ggtitle("Kendall correlations between university rankings",
          "CMDS correlation monoplot, second & third principal coordinates")

