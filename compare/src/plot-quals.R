# make supplemental figures for paper
options(stringsAsFactors=FALSE)
library(ggplot2)

#X11.options(antialias="subpixel", type="cairo")

args = commandArgs(TRUE)

df = read.delim(args[1])
out_eps_or_png = args[2]

df[df$method == "bis2", "method"] = "bismark-bt2"
df[df$method == "bis1", "method"] = "bismark-bt1"
df[df$method == "bwastrand", "method"] = "bwa-strand"
df[df$method == "bwa", "method"] = "bwameth"

df[grep("sim_R1", df$method, fixed=TRUE), "method"] = "bison"
df[grep("real_R1", df$method, fixed=TRUE), "method"] = "bison"
df$method = factor(df$method)#, levels=c("last", "bsmap", "gsnap", "bwa", "bismark", "bwa-strand", "bsmooth", "bison"))
df$method = factor(as.character(df$method), levels=rev(levels(df$method)))
df = df[df$qual > 0,]

df = df[order(df$qual),]

df = df[(df$on + df$off) != 0,]
df$on = df$on * 100
df$off = df$off * 100

require("grid")
print(dim(df))

SIM=any(grep("sim", args))
#         geom_point(color="grey60", size=ifelse(any(grep("sim", args)), 1.26, 0.2)) +

p = ggplot(df, aes(x=off, y=on, by=method)) +
         geom_line(aes(color=method), linetype="dashed", linewidth=ifelse(SIM, 0.1, 0.25) ) +
         geom_point(aes(color=method), size=ifelse(SIM, 1.0, 0.55)) +
         scale_shape(solid = FALSE) 
p = p + ylab("% Reads On Target")
p = p + xlab("% Reads Off Target")
p = p + scale_color_brewer(palette="Set1")
if(any(grep("sim", args))){
    p = p + xlim(xmin=0, xmax=1.5)
    p = p + ylim(ymin=65, ymax=92)
    p = p + geom_point(data=df[df$method %in% c('bismark-bt1', 'bismark-bt2', 'gsnap', "bsmap"),],
                       aes(x=off, y=on, by=method, color=method))
} else {
    # make the points larger, but don't try this for the bwa strand comparison
    if(!any(grep("strand", args))){
        p = p + geom_point(data=df[df$method %in% c('bismark-bt1', 'bismark-bt2', 'gsnap', "bsmap"),],
                       aes(x=off, y=on, by=method, color=method))
    }
    p = p + xlim(xmin=4, xmax=14)
    p = p + ylim(ymin=55, ymax=80)
}
#p = p + theme_bw()
p = p + theme(
             #legend.position = c(0.55, 0.25),
             #legend.position = c(ifelse(SIM, 0.69, 0.34), 0.25),
             #legend.position = ifelse(SIM, c(0.69, 0.25), c(0.29 ,0.75)),
             legend.position = c(0.71 ,0.25),
              legend.text=element_text(size=5, lineheight=4),
              axis.text=element_text(size=5),
              axis.title=element_text(size=7),
              legend.key.size=unit(5, "mm")

#              legend.key.height=5
              )
if(any(grep("strand", args))){
   p = p + theme(
             legend.position = c(ifelse(SIM, 0.69, 0.34), 0.25)
    )
}
p = p + guides(color=guide_legend(ncol=2, title=NULL))
#print(p)
ggsave(file=out_eps_or_png, units="cm", width=8.6, height=6.3,
    dpi=1200)

