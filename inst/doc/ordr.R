## ----knitr options, include=FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center", fig.width = 6, fig.height = 5
)

## ----setup--------------------------------------------------------------------
data(HairEyeColor)
library(MASS)
library(ordr)

## -----------------------------------------------------------------------------
print(HairEyeColor)
plot(HairEyeColor)

## -----------------------------------------------------------------------------
haireye <- apply(HairEyeColor, c(1L, 2L), sum)
haireye_ca <- corresp(haireye, nf = 3L)
print(haireye_ca)
# proportion of variance in each dimension
haireye_ca$cor^2 / sum(haireye_ca$cor^2)

## -----------------------------------------------------------------------------
# correspondence matrix (matrix of relative frequencies)
(haireye_p <- haireye / sum(haireye))
# row and column weights
(haireye_r <- rowSums(haireye) / sum(haireye))
(haireye_c <- colSums(haireye) / sum(haireye))
# matrix of standardized residuals
(haireye_s <-
    diag(1 / sqrt(haireye_r)) %*%
    (haireye_p - haireye_r %*% t(haireye_c)) %*%
    diag(1 / sqrt(haireye_c)))
# singular value decomposition
haireye_svd <- svd(haireye_s)
# row and column standard coordinates
diag(1 / sqrt(haireye_r)) %*% haireye_svd$u[, 1:3]
diag(1 / sqrt(haireye_c)) %*% haireye_svd$v[, 1:3]

## ---- fig.height=6------------------------------------------------------------
biplot(
  haireye_ca, type = "symmetric", cex = .8,
  main = "Correspondence analysis of subjects' hair & eye colors"
)

## -----------------------------------------------------------------------------
(haireye_ca_ord <- as_tbl_ord(haireye_ca))

## -----------------------------------------------------------------------------
get_conference(haireye_ca_ord)
confer_inertia(haireye_ca_ord, c(.25, .75))
confer_inertia(haireye_ca_ord, c(1, 1))
(haireye_ca_ord <- confer_inertia(haireye_ca_ord, "symmetric"))

## -----------------------------------------------------------------------------
glance(haireye_ca_ord)

## -----------------------------------------------------------------------------
augment_ord(haireye_ca_ord)

## ----tidy---------------------------------------------------------------------
tidy(haireye_ca_ord)

## ----scree plot---------------------------------------------------------------
ggplot(tidy(haireye_ca_ord), aes(x = name, y = inertia)) +
  geom_col() +
  labs(x = "Component", y = "Inertia") +
  ggtitle("Correspondence analysis of subjects' hair & eye colors",
          "Decomposition of inertia")

## ----fortify------------------------------------------------------------------
fortify(haireye_ca_ord)

## -----------------------------------------------------------------------------
haireye_ca_ord %>%
  augment_ord() %>%
  fortify() %>%
  transform(feature = ifelse(.matrix == "rows", "Hair", "Eye")) %>%
  ggbiplot(aes(color = feature, shape = feature, label = name), clip = "off") +
  theme_biplot() +
  geom_origin() +
  geom_rows_point() +
  geom_cols_point() +
  geom_rows_text(vjust = -1, hjust = 0, size = 3) +
  geom_cols_text(vjust = -1, hjust = 0, size = 3) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  scale_size_area() +
  ggtitle("Correspondence analysis of subjects' hair & eye colors",
          "Symmetric biplot")

## -----------------------------------------------------------------------------
sessioninfo::session_info()

