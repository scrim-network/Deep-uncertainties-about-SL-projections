# Analyse affect interpretation of likely range on projected probabilistic range

# clear environment
rm(list = ls())
graphics.off()

library(reshape2)
library(lhs)

library(ggplot2)
library(grid)
library(gridExtra)
library(rriskDistributions)

# Function to plot two plots along the same axis
source("ggplot_dual_axis.R")
# different strategies to estimate parameters from median and higher quantile
# and sample
source("sampling_strategies.R")

# read projected ice sheet ranges
df <- read.csv("../../data/Bamber2013.csv")
df <- df[which(df$Elicitation==2012 & df$Sheet=="WAIS"),-3:-4] # sample WAIS and elicitation 2012

# ================================================================================================
# integrate to cumulative SLR contribution, assuming linear increase of melt rate
# and zero [mm/yr] WAIS contribution around 1995
df[3:5] <- df[3:5] * (2090-1995)/(2100-1995)      # rate at 2090 assuming linear increase of rate
df[3:5] <- ((2090-1995) * df[3:5] / 2 ) / 1000    # factor 1000 to convert from [mm] to [m]

# Extract row-numbers with 'expert' estimates
experts <- which(df$agg==FALSE)                   

# Estimate median of expert quantile estimates
df <- rbind(df,
            data.frame(expert="Median",agg=TRUE,
                       p05=median(df$p05[experts]),
                       p50=median(df$p50[experts]),
                       p95=median(df$p95[experts])))

# Estimate max range
df <- rbind(df,
            data.frame(expert="MaxRange",agg=TRUE,
                       p05=min(df$p05[experts]),
                       p50=median(df$p50[experts]),
                       p95=max(df$p95[experts])))

pp <- (1:9999)/10000   # considered non-exceedance probabilities
pl = 0.05              # non-exceedance probability, lower quantile estimate
ph = 0.95              # non-exceedance probability, upper quantile estimate
mx = 3.3               # maximum ice sheet contribution WAIS

# Estimate cumulative distribution function fits
cdf   <- data.frame(expert=NULL,agg=NULL,pp=NULL,qq=NULL) # lognormal based
cdf.b <- data.frame(expert=NULL,agg=NULL,pp=NULL,qq=NULL) # beta      based

# only necessary for plotting
df.b  <- df

# fit cdfs to quantiles and use them to estimate quantiles of pp (non-exceedance probabilities)
for(i in 1:nrow(df)) {

  # expert and agg-expert quantiles
  ql  <- df$p05[i]
  q50 <- df$p50[i]
  qh  <- df$p95[i]
  
  # fit lognormal distribution
  if(ql==q50 || qh==q50) {
    qq <- rep(q50,length(pp))
    if(ql!=q50) {
      qq[which(pp<0.5)] <- sample.normal(q50, ql, pl, pp[which(pp<0.5)])
    }
    if(qh!=q50) {
      qq[which(pp>0.5)] <- sample.normal(q50, qh, ph, pp[which(pp>0.5)])
    }
} else {
    qq  <- sample.lnormal(ql=ql,q50=q50,qh=qh,ph=0.95,pl=0.05,pp=pp)
  }
  
  cdf <- rbind(cdf,
               data.frame(expert=df$expert[i],agg=df$agg[i],pp=pp,qq=qq) )
  
  # fit beta distribution
  qq    <- sample.beta(ql=ql,q50=q50,qh=qh,ph=0.95,pl=0.05,pp=pp,mn.l=-0.1,mx=mx)
  cdf.b <- rbind(cdf.b,
              data.frame(expert=df$expert[i],agg=df$agg[i],pp=pp,qq=qq) )
  
}

# row-number of individual expert based quantile estimates
exps <- which(!cdf$agg)

# add equal weight integrated estimate 'Integrated', i.e. combine all quantiles of all experts
# lnormal
#cdf <- rbind(cdf,
#             data.frame(expert="Integrated",agg=TRUE,pp=pp,
#                        qq=quantile(cdf$qq[exps],pp)))
#df <- rbind(df,
#            data.frame(expert="Integrated",agg=TRUE,
#                       p05=quantile(cdf$qq[exps],0.05),
#                       p50=quantile(cdf$qq[exps],0.50),
#                       p95=quantile(cdf$qq[exps],0.95)) )

# beta
#cdf.b <- rbind(cdf.b,
#               data.frame(expert="Integrated",agg=TRUE,pp=pp,
#                          qq=quantile(cdf.b$qq[exps],pp)))
#df.b <- rbind(df.b,
#              data.frame(expert="Integrated",agg=TRUE,
#                         p05=quantile(cdf.b$qq[exps],0.05),
#                         p50=quantile(cdf.b$qq[exps],0.50),
#                         p95=quantile(cdf.b$qq[exps],0.95)) )

# add uniform distribution (least informative) to cdf and cdf.b
cdf <- rbind(cdf,
             data.frame(expert="Uniform",agg=TRUE,pp=pp,
                        qq=mx*pp))
cdf.b <- rbind(cdf.b,
               data.frame(expert="Uniform",agg=TRUE,pp=pp,
                          qq=mx*pp))

df <- rbind(df,
            data.frame(expert="Uniform",agg=TRUE,
                       p05=0.05*mx,
                       p50=0.50*mx,
                       p95=0.95*mx))

df.b <- rbind(df.b,
              data.frame(expert="Uniform",agg=TRUE,
                         p05=0.05*mx,
                         p50=0.50*mx,
                         p95=0.95*mx))


## select probabilities for plotting
pp <- c(min(pp), (1:50)/100,  (510:900)/1000, (9100:9990)/10000)
pp.rows <- which(!is.na(match(cdf$pp, pp)))

cdf   <- cdf[pp.rows,]
cdf.b <- cdf.b[pp.rows,]

# ids of aggregated expert based quantile estimates (like exps)
exps <- which(!cdf$agg)
aggs <- which(cdf$agg)

# df.melt (like melt function), necessary for easy plotting
df.m <- rbind(data.frame(expert=df$expert,
                         agg=df$agg,
                         pp=0.05,
                         qq=df$p05),
              data.frame(expert=df$expert,
                         agg=df$agg,
                         pp=0.50,
                         qq=df$p50),
              data.frame(expert=df$expert,
                         agg=df$agg,
                         pp=0.95,
                         qq=df$p95))

df.m$expert <- levels(df.m$expert)[as.numeric(df.m$expert)]
df.m$expert[which(df.m$agg==FALSE)] <- "Expert"
df.m$expert <- factor(df.m$expert, levels=unique(df.m$expert))

# df.b.melt (like melt function)
df.b.m <- rbind(data.frame(expert=df$expert,
                           agg=df$agg,
                           pp=0.05,
                           qq=df$p05),
                data.frame(expert=df$expert,
                           agg=df$agg,
                           pp=0.50,
                           qq=df$p50),
                data.frame(expert=df$expert,
                           agg=df$agg,
                           pp=0.95,
                           qq=df$p95))

df.b.m$expert <- levels(df.m$expert)[as.numeric(df.m$expert)]
df.b.m$expert[which(df.m$agg==FALSE)] <- "Expert"
df.b.m$expert <- factor(df.m$expert, levels=unique(df.m$expert))

# ================================================================================================
# define colors, sizes and shapes, types used in plots
# ================================================================================================
colors <- c("Expert"     = "black",
            "EqualWts"   = "green",
            "PerfWts"    = "purple",
            "Median"     = "red",
            "MaxRange"   = "orange",
            "Integrated" = "black",
            "Uniform"    = "blue")

sizes <- c(2,5,0.1)
sizes <- c("Expert"     = sizes[2],
           "EqualWts"   = sizes[2],
           "PerfWts"    = sizes[2],
           "Median"     = sizes[2],
           "MaxRange"   = sizes[2],
           "Integrated" = sizes[3],
           "Uniform"    = sizes[3])

shapes <- c(4,16,NA)
shapes <- c("Expert"    = shapes[1],
           "EqualWts"   = shapes[2],
           "PerfWts"    = shapes[2],
           "Median"     = shapes[2],
           "MaxRange"   = shapes[2],
           "Integrated" = shapes[3],
           "Uniform"    = shapes[3])

types <- c("solid","dashed")
types <- c("Expert"      = types[1],
            "EqualWts"   = types[1],
            "PerfWts"    = types[1],
            "Median"     = types[1],
            "MaxRange"   = types[1],
            "Integrated" = types[1],
            "Uniform"    = types[2])

# ================================================================================================
## PLOTTING log-normal based quantile-functions ##
# ================================================================================================
# function to reverse log-scale axis)
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

label.positions <- c(1.0,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001)
breaks <- c(1,(9:1)/10,(9:1)/100,(9:1)/1000)
labels <- rep("",length(breaks))

label.positions <- 
  which(!is.na(match(breaks,label.positions)))
labels[label.positions] <- 
  breaks[label.positions]

p5 <- ggplot(data=cdf[aggs,], aes(x=1-pp, y=qq)) +
  scale_color_manual("",values=colors) +
  scale_size_manual("",values=sizes) +
  scale_shape_manual("",values=shapes) +
  scale_linetype_manual("",values=types) +
  geom_point(data=df.m, aes(color=expert, shape=expert, linetype=expert), size=3) +
  geom_line(data=cdf[exps,], aes(group=expert, color="Expert"), size=0.8) +
  geom_line(aes(color=expert, linetype=expert), size=1.5) +
  geom_point(data=df.m, aes(color=expert, shape=expert, size=expert)) +
  scale_x_continuous(trans=reverselog_trans(10),
                     breaks=breaks,
                     labels=labels) +
  scale_y_continuous(limits=c(NA,4)) +
  labs(x="Probability of exceedance [-]",
       y="West-Antarctic ice sheet contribution to sea level rise [m]\n(1995-2090)") +
 theme(legend.position=c(0.1,0.8),
        legend.key.height = unit(1,"lines"),
        legend.key=element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill=alpha('white', 0.0)), 
       panel.background = element_blank(),                                    # clear background color
       panel.grid.major = element_blank(),                                    # remove major grid
       panel.grid.minor = element_blank(),                                    # remove major grid
       panel.border = element_rect(colour = "black", fill=NA, size=1))        # add border around the plot

# Save quantile pdf
pdf("../../figures/quantile_functions_expert_judgments.pdf", paper="special", width=8, height=8/1.618)
print(p5)
dev.off()

# ================================================================================================
## PLOT SURVIVAL FUNCTIONS ##
# ================================================================================================
label.positions <- c(0.001,0.01,0.1,1)
breaks          <- c(0.001*(1:9), 0.01*(1:9), 0.1*(1:9), 1)
empty.labels <- rep("",length(breaks))

label.positions <- 
  which(!is.na(match(breaks,label.positions)))
labels <- empty.labels
labels[label.positions] <- 
  breaks[label.positions]

# basis plot
p.base <- ggplot(data=cdf[aggs,], aes(y=1-pp, x=qq)) +       # default aestetics
  geom_vline(xintercept=0:4, size=0.25, color='grey', linetype='dashed') +
  geom_hline(yintercept=c(0.001,0.01,0.1,1), size=0.25, color='grey', linetype='dashed') +
  scale_color_manual("",values=colors) +                     # scales
  scale_size_manual("",values=sizes) +
  scale_shape_manual("",values=shapes) +
  scale_linetype_manual("",values=types) +
  scale_x_continuous(limits=c(NA,4)) +                       # axes
  scale_y_log10(limits=c(0.001,1),
                breaks=breaks,
                minor_breaks=label.positions,
                labels=labels) +
  labs(x="West-Antarctic ice sheet contribution to sea level rise [m]\n(1995-2090)",                  # labels
       y="Probability of exceedance [-]") +
  theme(legend.position='right',                          # legend for 1  column journal
      legend.key.height = unit(1.5,"lines"),                 # legend
      legend.key=element_blank(),
      legend.text = element_text(size = 12),
      legend.background = element_rect(fill=alpha('white', 0.0)) )


# plot log_normal based estimates
p1 <- p.base +
  geom_point(data=df.m,      aes(color=expert, shape=expert, linetype=expert), size=2) +
  geom_line(data=cdf[exps,], aes(group=expert, color="Expert"), size=0.8) +
  geom_line(data=cdf[aggs,], aes(color=expert, linetype=expert), size=0.8) +
  geom_point(data=df.m,      aes(color=expert, shape=expert), size=2)

p2 <- p1 +
  labs(x=NULL) +
  theme(axis.title.x = element_text(margin = margin(t=0, b=0))) +
  scale_x_continuous(breaks = 0:4,         # axes
                     labels = rep("",5),
                     limits=c(NA,4))

g1 <- ggplot_dual_axis(p1, p2, which.axis = "x")
grid.draw(g1)

# plot beta based estimates
p3 <- p.base +
  geom_point(data=df.b.m,       aes(color=expert, shape=expert, linetype=expert), size=2) +
  geom_line( data=cdf.b[exps,], aes(group=expert, color="Expert"), size=0.8) +
  geom_line( data=cdf.b[aggs,], aes(color=expert, linetype=expert), size=0.8) +
  geom_point(data=df.b.m,       aes(color=expert, shape=expert), size=2) +
  theme(        panel.background = element_blank(),                                    # clear background color
                #panel.grid.major = element_blank(),                                    # remove major grid
                panel.grid.minor = element_blank(),                                    # remove major grid
                panel.border = element_rect(colour = "black", fill=NA, size=1))      # add border around the plot


p4 <- p3 +
  labs(x=NULL) +
  theme(axis.title.x = element_text(margin = margin(t=0, b=0))) +
  scale_x_continuous(breaks = 0:4,         # axes
                     labels = rep("",5),
                     limits=c(NA,4))

g2 <- ggplot_dual_axis(p3, p4, which.axis = "x")
grid.draw(g2)

# Save log normal survival function pdf 
pdf("../../figures/survival_functions_lnormal_expert_estimates.pdf", paper="special", width=5, height=6/1.618)
grid.draw(g1)
dev.off()

# for 2 column journal
#pdf("../../figures/survival_functions_beta_expert_estimates.pdf", paper="special", width=5, height=7/1.618)
#print(p3)#grid.draw(g2)
#dev.off()

# ================================================================================================
# Save as pdf (figure 3)
# ================================================================================================
# for 1 column journal
pdf("../../figures/survival_functions_beta_expert_estimates.pdf", paper="special", width=7.5, height=6/1.618)
print(p3)#grid.draw(g2)
dev.off()

postscript("../../figures/survival_functions_beta_expert_estimates.eps", paper="special", width=7.5, height=6/1.618)
print(p3)#grid.draw(g2)
dev.off()
# ================================================================================================
              
