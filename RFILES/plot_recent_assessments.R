# Change ranges for rcp8.5 (except rcp4.5 and period 1995-2090)

# clear environment
rm(list = ls())
graphics.off()

# set up component names
components <- c("GSL", "TE", "GSIC", "GIS", "AIS", "land")
fullnames  <- c("Global\nsea level", "Thermosteric\nexpansion", "Glaciers &\nsmall ice caps","Greenland\nice sheet","Antarctic\nice sheet","Water on/in\nsolid earth")

# read projected GSL and (component) ranges
df         <- read.csv("../../data/assessed_ranges.csv")
start.rate <- read.csv("../../data/rate0.csv")

t0 = 1995 # start year
tt = 2090 # sight year / time horizon
time.span  <- tt - t0

# ================================================================================================
# scale given ranges to match the time.span
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


# ================================================================================================
# reshape data for plotting
# ================================================================================================
# remove projections that are not plotted
df      <- df[which(df$Scenario!="A1B"),]        # exclude SRES A1B

# Add projection name and scenario information
df$Projection <- df$Assessment
levels(df$Projection)[levels(df$Projection)=="IPCC"] <- "(IPCC, 2013)"
levels(df$Projection)[levels(df$Projection)=="KNMI"] <- "(KNMI, 2014)"
levels(df$Projection)[levels(df$Projection)=="Kopp"] <- "(Kopp et al., 2014)"
levels(df$Projection)[levels(df$Projection)=="Grinsted"] <- "(Grinsted et al., 2015)"
levels(df$Projection)[levels(df$Projection)=="Perrette"] <- "(Perrette et al., 2013)"
#df$Projection <- paste(df$Assessment,"-",df$Scenario," (",df$Range,")",sep="")
df$Projection <- factor(paste(df$Scenario,"-", df$Range, df$Projection, sep=" "))

# Seperate dataset 66% projections for which also 90% is available (Kopp and Grinsted)
# and the other projections (not 66% or Perrette)
df66 <- df[which(df$Range=="66%" & df$Assessment!="Perrette"),]
dff  <- df[which(df$Range!="66%" | df$Assessment=="Perrette"),] 

# re-arrange factor order
dff$Assessment <- factor(dff$Assessment, levels = unique(dff$Assessment))
dff$Component  <- factor(dff$Component,  levels = unique(dff$Component))
dff$Projection <- factor(dff$Projection, levels = unique(dff$Projection))
df66$Projection <- factor(df66$Projection, levels = unique(df66$Projection))

# Determine plotting positions bars on the basis of Assessment and Component
rel.distance <- 7.5
positions    <- as.numeric(dff$Projection) - mean(as.numeric(dff$Projection))
xlocations   <- rel.distance * as.numeric(dff$Component)
dff$aux      <- positions + xlocations

df66$aux[which(df66$Assessment=="Kopp")]     <- dff$aux[which(dff$Assessment=="Kopp")]
df66$aux[which(df66$Assessment=="Grinsted")] <- dff$aux[which(dff$Assessment=="Grinsted")]

# ================================================================================================
# create a ggplot
# ================================================================================================
# set colors
projection.colors <- c("pink", "gold", "red",  "green2", "green4",  "magenta", "purple", "cyan2")
names(projection.colors) <- c(levels(dff$Projection),
                              levels(df66$Projection))[c(1:3, 7, 4, 8, 5:6)]

# create a ggplot
p <- ggplot(data=dff, aes(x=aux, y=p50, ymin=p05, ymax=p95)) +
  geom_linerange( size=2, aes(fill=Projection, color=Projection)) + 
  geom_linerange(data=df66, size=2, aes(fill=Projection, color=Projection)) + 
  geom_point(size=1.5) + 
  scale_color_manual(values=projection.colors, breaks=names(projection.colors)) +
    labs(y="Sea-level contribution [m]\n(1995-2090)", x=NULL, color=NULL) +
    theme(legend.position=c(0.6,0.82),                          # legend
        legend.key.height = unit(1,"lines"),                 # legend
        legend.key=element_blank(),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        axis.text.x = element_text(hjust = 1, angle=60, size=14),
        panel.background = element_blank(),                                    # clear background color
        panel.grid.major = element_blank(),                                    # remove major grid
        panel.grid.minor = element_blank(),                                    # remove major grid
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +      # add border around the plot
  scale_x_continuous(breaks=unique(xlocations), labels=fullnames)

p1 <-   p + theme(legend.position="right")

# ================================================================================================
# Save as pdf (figure 1)
# ================================================================================================
# prepared for presentation
pdf("../../figures/presentation_lsl_projections.pdf", paper="special", width=10, height=7/1.618); print(p1); dev.off()

# prepared for 2 column journal
#pdf("../../figures/lsl_projections.pdf", paper="special", width=5, height=10/1.618); print(p); dev.off()
#postscript("../../figures/lsl_projections.eps", paper="special", width=5, height=10/1.618);
#print(p); dev.off()

# prepared for 1 column journal
pdf("../../figures/lsl_projections.pdf", paper="special", width=8, height=7/1.618); print(p1); dev.off()
postscript("../../figures/lsl_projections.eps", paper="special", width=8, height=7/1.618);
print(p1); dev.off()
# ================================================================================================


