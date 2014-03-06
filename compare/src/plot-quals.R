# make supplemental figures for paper
options(stringsAsFactors=FALSE)
library(ggplot2)

#X11.options(antialias="subpixel", type="cairo")

args = commandArgs(TRUE)

df = read.delim(args[1])
out_eps_or_png = args[2]

df[df$method == "bis2", "method"] = "bismark-bt2"
df[df$method == "bis1", "method"] = "bismark-bt1"
df[df$method == "bwar", "method"] = "bwa-strand"
df[df$method == "bwa", "method"] = "bwameth"

df[grep("sim_R1", df$method, fixed=TRUE), "method"] = "bison"
df[grep("real_R1", df$method, fixed=TRUE), "method"] = "bison"
df$method = factor(df$method)#, levels=c("last", "bsmap", "gsnap", "bwa", "bismark", "bwa-strand", "bsmooth", "bison"))
df = df[df$qual > 0,]

df = df[order(df$qual),]

df = df[(df$on + df$off) != 0,]
df$on = df$on * 100
df$off = df$off * 100

require("grid")


p = ggplot(df, aes(x=off, y=on, by=method)) +
         geom_point(color="grey60", size=1.26) +
         geom_point(aes(color=method), size=1.0) +
         scale_shape(solid = FALSE) +
         geom_line(aes(color=method), linetype="dashed", linewidth=0.1, size=0.1) 
p = p + ylab("% Reads On Target")
p = p + xlab("% Reads Off Target")
if(any(grep("sim", args))){
    p = p + xlim(xmin=0, xmax=1.5)
    p = p + ylim(ymin=65, xmax=92)
}
#p = p + theme_bw()
p = p + theme(
             #legend.position = c(0.55, 0.25),
             legend.position = c(0.69, 0.25),
              legend.text=element_text(size=6, lineheight=5),
              axis.text=element_text(size=6),
              axis.title=element_text(size=8),
              legend.key.size=unit(6, "mm")

#              legend.key.height=5
              )

p = p + guides(color=guide_legend(ncol=2, title=NULL))
#print(p)
ggsave(file=out_eps_or_png, units="cm", width=8.6, height=6.3,
    dpi=1200)

