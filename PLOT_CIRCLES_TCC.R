library(circlize)
setwd("/home/vinicius/Documents/RESEARCH/LABISISMI/2021_GENOME_INSTABILITY/01_ASSEMBLY/OnlyCoordsWithCov")

# NC_002608.1 1 365425
# NC_002607.1 1 2014239
# NC_001869.1 1 191346

GetPlot = function(replicon, size, cols){
        # read coverage files
        SC = read.table(paste0("bc02_vs_OLDREF.bedgraph.", replicon))
        bc02 = read.table(paste0("bc05_vs_OLDREF.bedgraph.", replicon))
        bc03 = read.table(paste0("bc03_vs_OLDREF.bedgraph.", replicon))
        bc04 = read.table(paste0("bc06_vs_OLDREF.bedgraph.", replicon))
        
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


cols = c("#dc322eff", "#2e4273ff", "#1b9175ff", "#ed866fff", "#bebebeff")
sim = read.table("sim_vs_OLDREF.bedgraph.NC_001869.1")


colnames(sim) = c("chr", "start", "end", "value1")



circos.clear()


circos.track(ylim=c(0,1),panel.fun=function(x,y) {
        chr=CELL_META$sector.index
        xlim=CELL_META$xlim
        ylim=CELL_META$ylim
        circos.text(mean(xlim),mean(ylim),chr)
},bg.col=Dark2,bg.border=F,track.height=0.06)


circos.genomicTrack(data=df,panel.fun=function(region,value,...) {
        circos.genomicLines(region,value,type="l",col="grey50",lwd=0.6)
        circos.segments(x0=0,x1=max(chr_size),y0=100,y1=100,lwd=0.6,lty="11",col="grey90")
        circos.segments(x0=0,x1=max(chr_size),y0=300,y1=300,lwd=0.6,lty="11",col="grey90")
        #circos.segments(x0=0,x1=max(ref$V2),y0=500,y1=500,lwd=0.6,lty="11",col="grey90")
},track.height=0.08,bg.border=F)
circos.yaxis(at=c(100,300),labels.cex=0.25,lwd=0,tick.length=0,labels.col=col_text,col="#FFFFFF")
