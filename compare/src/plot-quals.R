options(stringsAsFactors=FALSE)
library(ggplot2)

df = read.delim('sim-qual-summary.txt')

df[df$method == "bis1", "method"] = "bismark"
df = df[df$qual > 0,]

#df$method = rep(c("bismark", "bsmap", "bwa", "gsnap", "last"), 256)

#df = df[(df$qual > 1) & (df$qual < 255),]

df = df[(df$on + df$off) != 0,]
df$on = df$on * 100
df$off = df$off * 100

require("grid")


p = ggplot(df, aes(x=off, y=on, by=method)) +
         geom_point(aes(color=method)) +
         geom_line(aes(color=method), linetype="longdash") 
        # geom_vline(xintercept=df[df$method == "bwa" & df$qual == 60, "off"], 
        #           linetype="dashed")
p = p + ylab("% Reads On Target")
p = p + xlab("% Reads Off Target")
p = p + theme(
             legend.position = c(0.55, 0.25),
              legend.text=element_text(size=5, lineheight=5),
              axis.text=element_text(size=5),
              axis.title=element_text(size=8),
              legend.key.size=unit(6, "mm")

#              legend.key.height=5
              )
p = p + labs(color="")
p = p + guides(color=guide_legend(ncol=2, title=NULL))
#print(p)
ggsave(file='qual-plot-real.eps', units="cm", width=8.6, height=6.3,
    dpi=600)

