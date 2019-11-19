#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Plots for LaTeX

library(copula)
library(ggplot2)
library(data.table)


# Copula plots ------------------------------------------------------------

copulas <- c("Frank", "Clayton", "Gumbel", "t")
u <- vector("list", length(copulas))
n <- 1000

set.seed(314)
for(i in seq_along(u)) {
  cop <- switch(copulas[i],
                "Frank" = frankCopula(5),
                "Clayton" = claytonCopula(2),
                "Gumbel" = gumbelCopula(2),
                "t" = tCopula(0.8, 2))
  u[[i]] <- data.table(rCopula(n, cop))
  u[[i]][, fam := copulas[i]]
}

u <- rbindlist(u)
setcolorder(u, "fam")
colnames(u)[2:3] <- c("u1", "u2")

u[, c("x1", "x2") := .(qnorm(u1), qnorm(u2))]

cop_fams <- ggplot(u, aes(x = x1, y = x2, group = fam)) +
  geom_point() + 
  #geom_density_2d() +
  facet_wrap(~fam) +
  theme_minimal() +
  xlim(-4, 4) +
  ylim(-4, 4)


ggsave("../tex/figures/cop_fams.pdf",
       cop_fams,
       device = "pdf",
       width = 10, height = 10, units = "cm",
       scale = 1.5)
