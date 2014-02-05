# make figure 1 from paper
options(stringsAsFactors=FALSE)
library(ggplot2)

df = read.delim('qual-summary.txt')

df[df$method == "bis1", "method"] = "bismark"
df$method = factor(df$method, levels=c("last", "bsmap", "gsnap", "bwa", "bismark"))
df = df[df$qual > 0,]

df = df[order(df$qual),]

df = df[(df$on + df$off) != 0,]
df$on = df$on * 100
df$off = df$off * 100

require("grid")


p = ggplot(df, aes(x=off, y=on, by=method)) +
         geom_point(aes(shape=method), size=1.4) +
         scale_shape(solid = FALSE)
         #geom_line(aes(linestyle=method), size=1.4) + scale_shape(solid=FALSE)
         #geom_line(aes(color=method), linetype="dotted") 
p = p + ylab("% Reads On Target")
p = p + xlab("% Reads Off Target")
p = p + theme_bw()
p = p + theme(
             #legend.position = c(0.55, 0.25),
             legend.position = c(0.75, 0.25),
              legend.text=element_text(size=6, lineheight=5),
              axis.text=element_text(size=6),
              axis.title=element_text(size=8),
              legend.key.size=unit(6, "mm")

#              legend.key.height=5
              )
#p = p + guides(color=guide_legend(ncol=2, title=NULL))
p = p + guides(shape=guide_legend(ncol=2, title=NULL))
#print(p)
ggsave(file='qual-plot-real.eps', units="cm", width=8.6, height=6.3,
    dpi=1200)

