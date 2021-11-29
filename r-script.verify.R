library(reshape)
library(ggplot2)
library(stringr)

risk.factors = read.table('risk.factors.txt', sep='\t', header=T, stringsAsFactors=F)
nature.2 = read.table('val.Nature.suppl.Table.2.txt', sep='\t', header=T, stringsAsFactors=F, comment='')  ## proteome

nature.2.sub = nature.2[which(nature.2$P.value > 0.05), c('Gene.Symbol', 'Ratio.24h', 'P.value.24h')]
colnames(nature.2.sub) = c('gene', 'ratio', 'p')
found = merge(nature.2.sub, risk.factors, by.x='gene', by.y='GeneID', all=F)
length(unique(found$gene))  ## 92 
found.up = found[which(found$p < 0.05 & found$ratio > 0), ]
length(unique(found.up$gene))  ## 8
found.down = found[which(found$p < 0.05 & found$ratio < 0), ]
length(unique(found.down$gene))  ## 8

found.up.pos = found.up[which(found.up$Coefficient > 0), ]
found.down.neg = found.down[which(found.down$Coefficient < 0), ]
length(unique(found.up.pos$gene))  ## 3
length(unique(found.down.neg$gene))  ## 6
2^(found.up.pos$ratio) - 1
1 - 2^(found.down.neg$ratio)

## fig. 5 ************
cns = colnames(nature.2)
c.idx = grep('Control.24h', cns)
v.idx = grep('Virus.24h', cns)
expr = c()
# gg = unique(found.up.pos$gene); cc = c('black', 'red')
gg = unique(found.down.neg$gene); cc = c('black', 'blue')
for(g in gg) {
	i = which(nature.2$Gene.Symbol == g)[1]
	temp.df = data.frame(ms=as.numeric(nature.2[i, c(c.idx, v.idx)]), group=c(rep('Control', length(c.idx)), rep('Virus', length(v.idx))), stringsAsFactors=F)
	temp.df$gene = g;
	expr = rbind(expr, temp.df)
}
expr$x = paste(expr$gene, expr$group, sep='\n')
expr.wk = expr
expr.wk[which(expr.wk$group == 'Control'), 'group'] = 'Ctrl'
expr.wk[which(expr.wk$group == 'Virus'), 'group'] = 'Infected'
expr.wk$x = paste(expr.wk$gene, expr.wk$group, sep='\n')
ggplot(expr.wk) + geom_boxplot(aes(x=x, y=ms, group=x, color=group)) + scale_y_continuous(trans='log10') + labs(x='', y='Protein Abundance') + theme_bw(base_size=18) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position='none') + scale_color_manual(values=cc) + geom_vline(xintercept=seq(2.5, length(gg)*2, 2), color='gray')

