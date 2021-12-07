ASM680v1library(circlize)

# NC_002608.1 1 365425
# NC_002607.1 1 2014239
# NC_001869.1 1 191346

# Function to plot with Circlize
GetPlot = function(replicon, size, cols){
        # read coverage files
        SC = read.table(paste0("NRC1_vs_ASM680v1.bedgraph.", replicon))
        bc02 = read.table(paste0("dura3_vs_ASM680v1.bedgraph.", replicon))
        bc03 = read.table(paste0("dsmap1_vs_ASM680v1.bedgraph.", replicon))
        bc04 = read.table(paste0("drnr_vs_ASM680v1.bedgraph.", replicon))
        
        # parse columns
        colnames(SC) = c("chr", "start", "end", "value1")
        colnames(bc02) = c("chr", "start", "end", "value1")
        colnames(bc03) = c("chr", "start", "end", "value1")
        colnames(bc04) = c("chr", "start", "end", "value1")
        
        chr_names = replicon
        chr_size = size
        
        # plot circlize
        col_text <- "grey40"
        circos.par(cell.padding = c(0.01, 0, 0.01, 0))
        circos.initialize(factors = chr_names,
                          xlim = matrix(c(rep(0,length(chr_names)), chr_size), ncol=2))
        
        # genomes
        circos.track(ylim=c(0,1),panel.fun=function(x,y) {
                chr=CELL_META$sector.index
                xlim=CELL_META$xlim
                ylim=CELL_META$ylim
                circos.text(mean(xlim),mean(ylim),chr)
        },bg.col="white",bg.border=F,track.height=0.021)
        
        skl = 10^3
        brk <- seq(0,400, 25)*skl
        circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
                circos.axis(h="top",major.at=brk,labels=round(brk/skl,1),
                            col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
        },bg.border=F)
        
        # cov
        circos.genomicDensity(SC, border = cols[1], col = cols[1], track.height = 0.1,bg.border=F)
        circos.genomicDensity(bc02, border = cols[2], col = cols[2], track.height = 0.1,bg.border=F)
        circos.genomicDensity(bc03, border = cols[3], col = cols[3], track.height = 0.1,bg.border=F)
        circos.genomicDensity(bc04,border = cols[4], col = cols[4], track.height = 0.1,bg.border=F)
        
}

GetPlot("NC_002608.1",
        365425,
        c("#2e4273ff", "#8491b4ff", "#1b9175ff", "#91d1c2ff"))

