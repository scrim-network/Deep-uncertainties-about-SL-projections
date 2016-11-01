# Analyse affect interpretation of likely range on projected probabilistic range

# clear environment
rm(list = ls())
graphics.off()

library(reshape2)
library(lhs)

library(ggplot2)
library(grid)
library(gridExtra)

#source("sampling_strategies.R")
source("~/Programs/R-scripts/set_theme_AB.R")
theme_set(theme_ali(base_size=25))
source("get_legend.R")

# ===============================================================
# IPCC likely ranges for period 1995-2090
A1B.likely = c(0.37,0.69) # A1B SRES
RCP.likely = c(0.45,0.82) # RCP8.5

v0 <- read.csv("../../data/rate0.csv")$GSL # start rate GSL-change [mm]

t0 = 1995 # start year
tt = 2090 # sight year / time horizon
time.span  <- tt - t0
# ===============================================================


# ==================================================================================================
# Horton et al. (2014) expert elicitation on change between 2000-2100 [cm]

t1 = 2000
t2 = 2100

# read expert opinions
horton <- read.csv("../../data/Horton_elicitation_RCP8.5_2100.csv") # [cm]
horton <- 10 * horton                                               # [mm]

# scale estimates for IPCC period
a      <- (horton + v0*(t1-t2) ) / ((t2 - t0)^2 - (t1 - t0)^2) # acceleration [mm/yr2]
horton <- (0 + v0 * time.span + a * time.span^2) / 1000        # change [m]

# share of experts that project a 17th percentile larger than the upper-bound of IPCC's likely range
# 27% of experts project a 83% probability that the sea-level rise will be larger than 
# 22 out of 81
P17s  <- horton$P17[which(!is.na(horton$P17))]
share <- length(which(P17s>RCP.likely[2])) / length(P17s)

# median
horton <- as.data.frame(t(apply(horton, 2, function(x) median(x, na.rm=T))))
# ==================================================================================================


# ==================================================================================================
# sem projections [meter]
sem <- read.csv("../../data/sem_ranges.csv")    # read data
sem <- sem[which(sem$Type=="SEM"),]             # and select SEM projections

# acceleration [mm/yr2]
a <- (1000*sem[,c("p05","p50","p95")] + v0*(sem$ref.year-sem$fut.year) ) /
  ((sem$fut.year - t0)^2 - (sem$ref.year - t0)^2)

# adjusted ranges [m]
sem[,c("p05","p50","p95")] <- (0 + v0 * time.span + a * time.span^2) / 1000

# ==================================================================================================


# ==================================================================================================
# combine sem and horton elicitation, and prepare for plotting
sem       <- rbind(sem,
                   data.frame(Assessment = "Horton et al. (2014)",
                              Type       = "EE95",
                              Scenario   = "RCP8.5",
                              p05        = horton$P05,
                              p50        = NA,
                              p95        = horton$P95,
                              ref.year   = 2000,
                              fut.year   = 2100) )

sem$id <- 1:nrow(sem)
A1B.ids <- which(sem$Scenario=="A1B")
RCP.ids <- which(sem$Scenario=="RCP8.5")
sem$id[RCP.ids] <- sem$id[RCP.ids] + 1

horton$id <- max(sem$id)
# ==================================================================================================


# ==================================================================================================
# set variables for plotting
split = max(A1B.ids) + 1

colors <- c("IPCC 'likely' range (RCP8.5)" = "gold",
            "IPCC 'likely' range (A1B)"    = "pink",
            "SEM"                          = "grey",
            "EE95"                         = "blue",
            "EE66"                         = "lightblue")
expand <- 0.6
A1B <- data.frame(x = c(min(A1B.ids)-expand,split),
                  y = A1B.likely)

RCP <- data.frame(x = c(split,max(sem$id[RCP.ids])+expand),
                  y = RCP.likely)

# plotting
p2 <- ggplot(data=sem, aes(x=id)) +
  geom_ribbon(data=A1B, aes(x=x, ymin=A1B$y[1], ymax=A1B$y[2], fill="IPCC 'likely' range (A1B)")) +
  geom_line(data=A1B, aes(x=x, y=A1B$y[1]), size=0.8, color="black") +
  geom_line(data=A1B, aes(x=x, y=A1B$y[2]), size=0.8, color="black") +
  geom_ribbon(data=RCP, aes(x=x, ymin=RCP$y[1], ymax=RCP$y[2], fill="IPCC 'likely' range (RCP8.5)")) +
  geom_line(data=RCP, aes(x=x, y=RCP$y[1]), size=0.8, color="black") +
  geom_line(data=RCP, aes(x=x, y=RCP$y[2]), size=0.8, color="black") +
  geom_vline(xintercept=split, color="black", size=1.2) +
  geom_linerange(aes(y=p50,       ymin=p05, ymax=p95, color=Type), size=6) +
  geom_linerange(data=horton, aes(ymin=P17, ymax=P83, color="EE66"), size=6) +
  scale_fill_manual("",values=colors) +
  scale_color_manual("",values=colors) +
  geom_point(aes(y=p50), color="black", size=4) +
  annotate(
    geom = "text", x = 1.5, y = 1.8, size = 8,
    label = "SRES A1B", hjust=1) +
  annotate(
    geom = "text", x = 10.5, y = 1.8, size = 8,
    label = "RCP8.5", hjust=1) +
  annotate(
    geom = "text", x = max(sem$id), y = -0.00, size = 4.5, color="blue",
    label = "(Expert elicited)", hjust=0) +
  expand_limits(y=0) +
  scale_y_continuous(limits=c(0,1.8)) +
  coord_flip() +
  scale_x_reverse(breaks=sem$id,
                  labels=sem$Assessment,
                  expand=c(0,0)) +
  labs(x = NULL, y = "Sea-level change [m]\n(1995 and 2090)") +
  theme(legend.position="none",
        panel.background = element_rect(fill = NULL, colour = "black", size = 2, linetype = 1),
        panel.ontop = TRUE,
        axis.title.y=element_text(angle=90, margin = margin(r=10)),
        axis.title.x=element_text(margin = margin(t=10, b=10)) )
    



pdf("../../figures/sem_projections.pdf", paper="special", width=10, height=15/1.618)
print(p2)
dev.off()

