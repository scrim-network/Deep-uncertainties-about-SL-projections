# Analyse affect interpretation of likely range on projected probabilistic range

# clear environment
rm(list = ls())
graphics.off()

library(reshape2)
library(lhs)

library(ggplot2)
library(grid)
library(gridExtra)

# different strategies to estimate parameters from median and higher quantile
# and sample
source("sampling_strategies.R")

# set up component names
components <- c("GSL", "TE", "GSIC", "GIS", "AIS", "land")

# read projected GSL and (component) ranges
df         <- read.csv("../../data/assessed_ranges.csv")
start.rate <- read.csv("../../data/rate0.csv")

t0 = 1995 # start year
tt = 2090 # sight year / time horizon
time.span  <- tt - t0

# ================================================================================================
# scale given ranges to match the time.span
# change ranges for rcp8.5 (except rcp4.5 and period 1995-2090)
# ================================================================================================
# approximate start rates [mm]
v0 <- df[,c("p05","p50","p95")]
for(comp in names(start.rate)) {
  v0[which(df$Component==comp),] <- start.rate[comp]
}

# acceleration [mm/yr2]
a <- (1000*df[,c("p05","p50","p95")] + v0*(df$ref.year-df$fut.year) ) /
  ((df$fut.year - t0)^2 - (df$ref.year - t0)^2)

# adjusted ranges [m]
df[,c("p05","p50","p95")] <- (0 + v0 * time.span + a * time.span^2) / 1000
# ================================================================================================

# Extract the 90% ranges from Kopp and KNMI
kopp90 <- as.numeric(df[which(df$Assessment=="Kopp" &
                                df$Range     =="90%"  &
                                df$Component =="GSL"), c("p05","p95")] )
knmi90 <- as.numeric(df[which(df$Assessment=="KNMI" &
                                df$Range     =="90%"  &
                                df$Component =="GSL"), c("p05","p95")] )

# Extract get process-based ranges as communicated by IPCC for RCP8.5 1990
GSL  <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="GSL"), c("p05","p50","p95")])
TE   <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="TE"),  c("p05","p50","p95")])
GSIC <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="GSIC"),c("p05","p50","p95")])
GIS  <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="GIS"),c("p05","p50","p95")])
AIS  <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="AIS"),c("p05","p50","p95")])
land <- as.numeric(df[which(df$Assessment=="IPCC" & df$Scenario=="RCP8.5" & df$Component=="land"),c("p05","p50","p95")])

# ================================================================================================
# estimate the 66, 90, and 99% range of GSL given different interpretations of the 'likely' range
# using sampling techniques
# ================================================================================================
n <- 25000
ranges <- 50:99/100
qq     <- c(0.005,0.05,0.13,0.50,0.83,0.95,0.995)

deep.unc <- as.data.frame(matrix(NA,nrow=length(ranges),ncol=length(qq)))
names(deep.unc) <- c("q005","q05","q17","q50","q83","q95","q995")

deep.unc <- cbind(data.frame(Range = ranges*100), deep.unc)

for(i in 1:length(ranges)) {
  ph <- 0.5 + ranges[i]/2
  pl <- 0.5 - ranges[i]/2
  
  df <- as.data.frame(randomLHS(n, 5))
  names(df) <- c("TE","GSIC","GIS","AIS","land")

  df$TE   <- sample.normal(q50=TE[2],qh=TE[3],ph=ph,pp=df$TE)               # TE   sampling from normal
  df$GSIC <- sample.normal(q50=GSIC[2],qh=GSIC[3],ph=ph,pp=df$GSIC)         # GSIC sampling from normal
  df$GIS  <- sample.lnormal(ql=GIS[1],q50=GIS[2],qh=GIS[3],ph=ph,pl=pl,pp=df$GIS) # GIS  sampling from lnormal
  df$AIS  <- sample.lnormal(ql=AIS[1],q50=AIS[2],qh=AIS[3],ph=ph,pl=pl,pp=df$AIS) # AIS  sampling from lnormal
  df$land <- sample.normal(q50=land[2],qh=land[3],ph=ph,pp=df$land)         # land sampling from normal

  deep.unc[i,-1] <- c(quantile(apply(df        , 1, sum), qq))
}

# ================================================================================================
# reshape data for plotting
# ================================================================================================
melt.deep <- rbind(
  data.frame(
    IPCC       = deep.unc$Range,
    Projection = "99%",
    central    = deep.unc$q50,
    lower      = deep.unc$q005,
    upper      = deep.unc$q995),
  data.frame(
    IPCC       = deep.unc$Range,
    Projection = "90%",
    central    = deep.unc$q50,
    lower      = deep.unc$q05,
    upper      = deep.unc$q95),
  data.frame(
    IPCC       = deep.unc$Range,
    Projection = "66%",
    central    = deep.unc$q50,
    lower      = deep.unc$q17,
    upper      = deep.unc$q83)
)

# ================================================================================================
# create a ggplot
# ================================================================================================
# set colors and color labels
clabels <-  c("Central estimate"       = "black",
              "66% range"              = "blue",
              "90% range"              = "deepskyblue",
              "99% range"              = "lightblue1",
              "IPCC 'likely' range"    = "gold",
              "90% (Kopp et al. 2014)" = "green4",
              "90% (KNMI'14)"          = "red",
              "Semi-empirical 90%"     = "grey",
              "High-end expert\n95% ice sheet contribution" = "orange",
              "Deep uncertainty due to\n hydrofracturing and cliff failure" = "grey")


colors <- c("Central estimate"    = "black",
            "66%"                 = "blue",
            "90%"                 = "deepskyblue",
            "99%"                 = "lightblue1",
            "IPCC 'likely' range" = "gold",
            "Kopp's 90% range"    = "green4",
            "KNMI's 90% range"    = "red",
            "Semi-empirical 90%"  = "grey",
            "High-end expert\n95% ice sheet contribution" = "orange",
            "Deep uncertainty due to\n hydrofracturing and cliff failure" = "grey"
            )

xbreaks <- c(50,60,66,70,80,90,100)
xlabels <- c("50","60","66\n(Kopp et al. 2014)","70","80","90\n(KNMI'14)","100")

# line thickness in y-units
yth = 0.015

# create a ggplot
p1 <- 
  ggplot(melt.deep[which(melt.deep$IPCC>=66),], aes(x=IPCC)) +
  scale_fill_manual("",values=colors) +
  scale_color_manual("", values=colors, breaks=names(colors), labels=names(clabels)) +

  geom_line(aes(y=(0.45+0.82)/2, color="66%"), size=0.8) +
  geom_line(aes(y=(0.45+0.82)/2, color="90%"), size=0.8) +
  geom_line(aes(y=(0.45+0.82)/2, color="99%"), size=0.8) +
  geom_ribbon(
    data=melt.deep[which(melt.deep$IPCC>=66),],
    aes(ymin=lower, ymax=upper, fill=Projection), alpha=1) +
  geom_line(
    data=melt.deep[which(melt.deep$IPCC>=66),],
    aes(y=central, color="Central estimate"), size=1) +
  geom_vline(xintercept=66, linetype="dashed") +
  geom_segment(aes(x=66,xend=66,y=kopp90[1],yend=kopp90[2], color="Kopp's 90% range"), size=2) +
  geom_segment(aes(x=66,xend=66,y=kopp90[2]-yth,yend=kopp90[2]), color=colors["Kopp's 90% range"], size=6) +
  geom_segment(aes(x=66,xend=66,y=kopp90[1],yend=kopp90[1]+yth), color=colors["Kopp's 90% range"], size=6) +
  geom_line(data=data.frame(x=c(90,90), y=knmi90), aes(x=x, y=y, color="KNMI's 90% range"), size=2) +
  geom_segment(aes(x=90,xend=90,y=knmi90[2]-yth,yend=knmi90[2]), color=colors["KNMI's 90% range"], size=6) +
  geom_segment(aes(x=90,xend=90,y=knmi90[1],yend=knmi90[1]+yth),color=colors["KNMI's 90% range"], size=6) +
  labs(y="Global sea level change [m]\n(1995-2090)", x="Probabilistic interpretation of\nIPCC 'likely' ranges [%]") +
  expand_limits(y=-0.1) +
  scale_y_continuous(limits=c(0,1.2)) +
  scale_x_continuous(limits=c(60,100), breaks=xbreaks, labels=xlabels) +
  guides(fill=FALSE) +
  theme(legend.position="right",                          # legend
        legend.key.height = unit(1,"lines"),                 # legend
        legend.key=element_blank(),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        panel.background = element_blank(),                                    # clear background color
        panel.grid.major = element_blank(),                                    # remove major grid
        panel.grid.minor = element_blank(),                                    # remove major grid
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +      # add border around the plot
  geom_segment(aes(x=66, xend=99, y=0.125, yend=0.125), size=1, arrow=arrow(angle = 15, ends = "both", type = "closed")) +
  geom_segment(aes(x=62, xend=62, y=0.45, yend=0.82), size=0.5, arrow=arrow(angle = 15, ends = "both", type = "closed"), color="gold") +
  geom_segment(aes(x=62, xend=62, y=0.5, yend=0.75), size=1,  color="gold") +
  geom_text( aes(x=82.5, y= 0.125, label="'likely' range spans\n66% probability or more")) +
  geom_label(aes(x=82.5, y= 0.125, label="'likely' range spans\n66% probability or more")) +
  geom_line(aes(y=0.45, color="IPCC 'likely' range"), size=1) +
  geom_line(aes(y=0.82, color="IPCC 'likely' range"), size=1)
  
#pdf("../../figures/likely_interpretation.pdf", paper="special", width=5, height=7.5/1.618)
#print( p1 )
#dev.off()

#postscript("../../figures/likely_interpretation.eps", paper="special", width=5, height=7.5/1.618)
#print(p1); dev.off()

# ================================================================================================
# Save as pdf (figure 2)
# ================================================================================================
pdf("../../figures/likely_interpretation.pdf", paper="special", width=7.5, height=6/1.618)
print( p1 )
dev.off()

postscript("../../figures/likely_interpretation.eps", paper="special", width=7.5, height=6/1.618)
print(p1); dev.off()
# ================================================================================================


