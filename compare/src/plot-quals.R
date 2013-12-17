library(ggplot2)

df = read.delim('qual-summary.txt')

#df$method = rep(c("bismark", "bsmap", "bwa", "gsnap", "last"), 256)

#df = df[(df$qual > 1) & (df$qual < 255),]

df = df[(df$on + df$off) != 0,]
df$on = df$on * 100
df$off = df$off * 100



p = ggplot(df, aes(x=off, y=on, by=method)) +
         geom_point(aes(color=method)) +
         geom_line(aes(color=method), linetype="longdash") +
         geom_vline(xintercept=df[df$method == "bwa" & df$qual == 60, "off"], 
                    linetype="dashed")
p = p + ylab("% Reads On Target")
p = p + xlab("% Reads Off Target")
#print(p)
ggsave('qual-plot.png')

