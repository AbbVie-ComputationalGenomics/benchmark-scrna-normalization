#library(tidyr)
library(ggplot2)
library(argparse)

cell_filter = function(sce, colstr, up_thresh_mad = TRUE, up_mad = 5, 
                       down_thresh_mad = TRUE, down_mad = 5,
                       up_fixed = NA, down_fixed = NA){
  # for filtering on total counts, number of genes, etc.
  df = data.frame(colData(sce))
  metric = df[,colstr]
  if(up_thresh_mad){
    up_thresh = median(metric, na.rm=T) + up_mad*mad(metric, na.rm=T)
    uptitle = paste0("> ",up_mad,"*mad")
  }else{
    up_thresh = up_fixed
    uptitle = paste0("> ",up_fixed)
  }
  if(down_thresh_mad){
    down_thresh = median(metric, na.rm=T) - down_mad*mad(metric, na.rm=T)
    downtitle = paste0("< ",down_mad,"*mad")
  }else{
    down_thresh = down_fixed
    downtitle = paste0("< ",down_fixed)
  }
  
  cellsbelow = length(which(metric<down_thresh))
  cellsabove = length(which(metric>up_thresh))
  
  cells_to_filter = rownames(df)[which(metric<down_thresh | metric>up_thresh)]
  nacells = rownames(df)[which(is.na(metric))]
  cells_to_filter = unique(c(cells_to_filter, nacells))
  numfilt = length(cells_to_filter)
  plot_title = paste0("Filtering cells with ",colstr," ",downtitle," or ",uptitle)
  subtitle = paste0(numfilt," (",round(100*numfilt/nrow(df), digits=2),"%) cells filtered")
  
  max = max(density(metric, na.rm=T)$y)
  ythresh = 0.75*max
  ycellnum = 0.25*max
  minm = min(metric, na.rm=T)
  maxm = max(metric, na.rm=T)
  
  p = ggplot(df, aes_string(colstr)) + geom_density() +
    geom_vline(xintercept = down_thresh, color="darkred") + 
    geom_vline(xintercept = up_thresh, color="darkred")+
    geom_label(data = data.frame(x=c(down_thresh, up_thresh), y=c(ythresh, ythresh), 
                                 label=round(c(down_thresh, up_thresh), digits=2)),
               aes(x=x,y=y,label=label), 
               inherit.aes=F, color="darkred") +
    annotate("text",x=(down_thresh-minm)/2, y=ycellnum, label=cellsbelow)+
    annotate("text",x=down_thresh+(up_thresh-down_thresh)/2, y=ycellnum, label=nrow(df)-numfilt)+
    annotate("text",x=up_thresh+(max(c(maxm,up_thresh))-up_thresh)/2, y=ycellnum,label=cellsabove)+
    annotate("text",x=0.05*maxm, y=max, label=paste0("NA: ",length(nacells))) +
    labs(title = plot_title, subtitle = subtitle) +
    theme(plot.subtitle = element_text(hjust = 0.5))
  
  summary = c(filter=colstr,up_thresh=up_thresh, down_thresh = down_thresh, 
              n_orig = nrow(df), n_failed = numfilt)
  return(list(cells_to_filter = cells_to_filter, plot= p, summary = summary))
}


gene_filter = function(sce, pct_cells= 5, min_cells = NA, colname="n_cells_by_counts"){
  n_cells = rowData(sce)[,colname]
  if(is.na(min_cells)){
    min_cells = floor(0.01*pct_cells*ncol(sce))
    plot_title = paste0("Filtering genes detected in less than ",
                        pct_cells, "% of cells")
  }else{
    plot_title = paste0("Filtering genes detected in less than ",
                        min_cells, " cells")
  }
  genes_keep = rownames(sce)[which(n_cells > min_cells)]
  n_genes_keep = length(genes_keep)
  n_failed = nrow(sce)-n_genes_keep
  summary = c(filter="genes_min_cells", 
              up_thresh = NA, 
              down_thresh = min_cells, 
              n_orig = nrow(sce), 
              n_failed = n_failed)
  max = max(density(n_cells, na.rm=T)$y)
  subtitle = paste0(n_failed," (",round(100*n_failed/nrow(sce), digits=2),
                    "%) genes filtered")
  
  p = ggplot(data.frame(rowData(sce)), aes_string(colname)) + 
    geom_density() + 
    geom_vline(xintercept = min_cells, color="darkred")+
    geom_label(data = data.frame(x=min_cells, y=0.75*max, 
                                 label=min_cells),
               aes(x=x,y=y,label=label), 
               inherit.aes=F, color="darkred") +
    annotate("text",x=min_cells/2, y=0.25*max, label=n_failed)+
    annotate("text",x=min_cells+(max(n_cells)-min_cells)/2, y=0.25*max, label=length(genes_keep))+
    labs(title = plot_title, subtitle = subtitle) + 
    theme(plot.subtitle = element_text(hjust=0.5))
  
  return(list(genes_to_keep = genes_keep, plot = p, summary = summary))
}

plot_img <- function ( Filtering, imgname ){
  imgpath <- paste0( 'filtering/', imgname )
  jpeg( imgpath, width=1000, height=800, res=200 )
  print (Filtering$plot)
  dev.off()
}


main <- function( dfile ){
  parser <- ArgumentParser( description='in house filtering function')
  parser$add_argument('--dfile', dest='dfile')
  parser$add_argument('--totalcount_up', dest='totalcount_up', type='integer')
  parser$add_argument('--totalcount_down', dest='totalcount_down', type='integer')
  parser$add_argument('--genenum_up', dest='genenum_up', type='integer')
  parser$add_argument('--genenum_down', dest='genenum_down', type='integer')
  parser$add_argument('--mitoperc_up', dest='mitoperc_up', type='integer')
  parser$add_argument('--mitoperc_down', dest='mitoperc_down', type='integer')
  parser$add_argument('--mincell_gene', dest='mincell_gene', type='integer')
  
  args <- parser$parse_args()
  
  sce = readRDS( args$dfile )
  #Filtering
  ## log10 total counts
  sce$log10sum = log10(sce$sum + 1)
  cf1 = cell_filter(sce, "log10sum", up_mad = args$totalcount_up, down_mad = args$totalcount_down)
  plot_img( cf1, 'filtering_log10sum.jpg' ) 
  
  ## log10 number of genes
  sce$log10detected = log10(sce$detected+1)
  cf2 = cell_filter(sce, "log10detected", up_mad = args$genenum_up, down_mad = args$genenum_down)
  plot_img( cf2, 'filtering_log10detected.jpg' )
  
  ## %mito
  cf3 = cell_filter(sce, "subsets_mito_percent", up_thresh_mad = FALSE, down_thresh_mad = FALSE, up_fixed = args$mitoperc_up, down_fixed= args$mitoperc_down)
  plot_img( cf3, 'filtering_mito_perc.jpg' )
  
  ## gene filter
  gf = gene_filter(sce, min_cells = args$mincell_gene, colname = "detected")
  plot_img( gf, 'filtering_genes.jpg' )
  
  #### add all filterings
  filtcells = unique(c(cf1$ccells_to_filter, cf2$cells_to_filter, cf3$cells_to_filter))
  print ( length(filtcells) )
  keepcells = colnames(sce)[which(!colnames(sce) %in% filtcells)]
  print (length(gf$genes_to_keep) )  ##TODO colnames(sce)
  
  sce_filt = sce[gf$genes_to_keep, keepcells]
  saveRDS(sce_filt, file="sce_filtered.RDS")
}

dir.create( "filtering" )
main( dfile )

