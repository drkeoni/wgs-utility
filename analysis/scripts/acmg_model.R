require(ggplot2)
require(plyr)

#
# used in cleaning up raw Excel spreadsheet (contains mix of string and
# floating point percentages)
#
convert.string.pct <- function(s) {
  i <- grep('%',s)
  f <- as.double(sub('%','',s))
  f[i] <- f[i] / 100.0
  f
}

#
# csv.file is something like 'exome_acmg56_model.cleaned.csv'
#
load.model <- function(csv.file) {
  d <- read.csv(csv.file)
  d2 <- ddply(d,c('condition'),function(x) {data.frame(gene.count=nrow(x))} )
  d <- merge(d,d2)
  d$per.gene.carrier.freq <- d$carrier.freq / d$gene.count
  d
}

shorten <- function(s,l) {
  sapply( as.character(s), function(a) {
    if (nchar(a)>l) {
      sprintf('%s...%s',substr(a,1,(l-3)/2),substr(a,nchar(a)-(l-3)/2+1,nchar(a)))
    } else {
      a
    }
  } )
}

#
# d is load.model
#
calculate.sensitivity <- function( d, snv.sens, indel.sens, cnv.sens ) {
  d$snv.sens <- snv.sens * d$adj.coverage
  d$indel.sens <- indel.sens * d$adj.coverage
  d$cnv.sens <- cnv.sens * d$adj.coverage
  d$weighted.sens <- d$snv.burden.est * d$snv.sens + 
                     d$indel.burden.est * d$indel.sens +
                     d$cnv.burden.est * d$cnv.sens
  d$prob.fn <- d$per.gene.carrier.freq * (1-d$weighted.sens)
  d
}

#
# d is load.model
#
plot.condition.carrier.freqs <- function(d) {
  d <- d[order(d$carrier.freq),]
  d$short.cond <- factor(shorten(d$condition,37),ordered=T,levels=shorten(d$condition,37))
  p <- qplot( short.cond, carrier.freq, data=d )
  p <- p + scale_x_discrete(name='Condition')
  p <- p + scale_y_log10(name='Carrier Frequency')
  p <- p + coord_flip()
  p
}
