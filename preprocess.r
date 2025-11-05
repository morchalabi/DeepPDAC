# ============================================================
# Author: Dr. Mori (Morteza) Chalabi
# Created: October 2025
# Description: Deep learning model trained on PDAC TME (LIV-MR signature)
# ============================================================

require(Seurat)
require(DoubletFinder)
require(ggplot2)
require(infotheo)
require(viridis)
require(Matrix)
set.seed(42)

# Function for identifying doublets ####
# 1: it creates n (pN) artificial doublets from processed Seurat object
# 2: it merges those doublets with original data and runs Seurat again
# 3: it finds an optimal k (pK) nearest neighbors of each original cell; it first tries several pN-PK choices
# 4: it computes proportion of simulated doublets (pANN: artificial nn) in neighborhood of each cell
# 5: plots histogram of pANN values
# 6: it selects a cutoff on histogram considering expected number of doublets (nExp) given by user
# 7: original cells with cutoff <= pANN are reported doublets
# Att.: nExp should be an estimate of doublets (Poisson doublet formation rate) adjusted for homotypic doublets
doubletIndentify = function(s_)
{
  # pK Identification (ground-truth ####
  
  sweep.res.list = paramSweep(seu = s_, PCs = 1:30, sct = FALSE, num.cores = 3)     # steps 1-3: trying different pN and pK
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)                          # computes BCmvn (Mean-variance normalized bimodality coefficient) needed to find optimal pK
  bcmvn_ = find.pK(sweep.stats)                                                     # step 3: finds optimal pK
  pk_ = as.numeric(as.character(bcmvn_$pK[which.max(bcmvn_$BCmetric)]))
  
  # homotypic Doublet Proportion Estimate ####
  
  nExp_poi = round(0.09*nrow(s_@meta.data))                         # assuming 9% doublet formation rate (from Poisson loading) including hetrotypic and homotypic
  homotypic.prop = modelHomotypic(s_@meta.data$seurat_clusters)     # ex: annotations = s_@meta.data$ClusteringResults
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))                 # without adjusting nExp_poi for homotypic doublets, all transitional cells are wrongly reported doublets leading to overestimation 
  
  # identifying doublets ####
  
  s_ = doubletFinder(s_, PCs = 1:30, pN = 0.25, pK = pk_, nExp = nExp_poi.adj, sct = FALSE)
  colnames(s_@meta.data)[grepl(pattern = 'DF.classification', x = colnames(s_@meta.data))] = 'doublet'      # renaming the doublet classification column
  
  return(s_$doublet)
}

# Function for identifying malignant cells ####
malignant_cells_identifier = function(s_, genes_)
{
  signature_len = length(genes_)                                    # number of genes in the malignancy signature
  genes_ = intersect(genes_, rownames(s_))                          # keeping only signature genes present in the input sample
  if(length(genes_) < 0.8*signature_len) warning('Less than 80% of malignant signature genes are present in this sample\n')
  
  avg_exprs = colMeans(s_[['RNA']]$data[genes_,,drop = F])          # average expression of signature genes in each cell
  
  min_avg_expr = quantile(avg_exprs[0 < avg_exprs], probs = .2)     # minimum average expression of signature genes in a cell to be considered malignant
  malignant_cells = colnames(s_)[min_avg_expr <= avg_exprs]         # supposedly malignant cells
  
  return(malignant_cells)
}

# Function for computing feature vector
f = function(...)
{
  # this function has been made hidden until publicaiton
  return(f_vect)
}

plotting = function(s_, cells_, sample_)
{
  p_ = list()
  
  ## plotting QC readouts ####
  
  # violin plot of gene counts colored for CD45 (immune marker)
  # immune cells (CD45+) are quite stable in terms of gene counts (200 <= nFeature_RNA <= 2500)
  # cancer-immune doublets will be CD45+ with very high (amplification) or low (loss) gene counts
  
  p_[[1]] = VlnPlot(s_, features = "nFeature_RNA", log = T, layer = 'counts', group.by = 'orig.ident')+
            theme(text = element_text(face = 'bold', size = 25), axis.text.y = element_text(size = 30), axis.text.x = element_blank() ,legend.position = 'none')+
            labs(x = 'Cells', y = 'Gene counts', title = '# Unique genes')+
            geom_hline(yintercept = c(200,5e3,10e3,15e3), color = 'black', alpha = 0.5, linetype = 'dashed')
  
  indx_ = discretize(X = s_[['RNA']]$counts['PTPRC',], nbins = 10)$X
  cols_ = colorRampPalette(colors = c('blue1','red4'))(10)
  cols_ = cols_[indx_]
  alphas_ = seq(0.3,1,length.out = 10)
  alphas_ = alphas_[indx_]
  sizes_ = seq(0.5,1.1,length.out = 10)
  sizes_ = sizes_[indx_]
  p_[[1]]$layers[[1]]$aes_params$fill = 'grey85'
  p_[[1]]$layers[[1]]$aes_params$colour = 'grey70'
  p_[[1]]$layers[[2]]$aes_params$colour = cols_
  p_[[1]]$layers[[2]]$aes_params$alpha = alphas_
  p_[[1]]$layers[[2]]$aes_params$size = sizes_
  
  # violin plot of percent of MT gene counts
  
  p_[[2]] = VlnPlot(s_, features = "percent.mt", ncol = 1, log = F, layer = 'counts', group.by = 'orig.ident')+
            theme(text = element_text(face = 'bold', size = 25), axis.text.y = element_text(size = 30), axis.text.x = element_blank() ,legend.position = 'none')+
            labs(x = 'Cells', y = 'Gene counts', title = 'Mitochondiral percentage')+
            geom_hline(yintercept = c(7, 10, 20), color = 'black', alpha = 0.5, linetype = 'dashed')
  
  # scatter plot of gene count vs MT percent
  
  p_[[3]] = ggplot(data = s_@meta.data, aes(y = nFeature_RNA, x = percent.mt, color = nFeature_RNA))+
            theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),
                  text = element_text(face = 'bold', size = 20),
                  axis.text.y = element_text(size = 15),axis.line = element_line(color = 'black'))+
            labs(x = 'Percent MT', y = 'Gene counts', title = NULL)+
            geom_point(size = 0.4, show.legend = F)+
            coord_transform(y = 'sqrt', x = 'sqrt', xlim = range(s_$percent.mt), ylim = range(s_$nFeature_RNA))+
            scale_color_viridis(option = 'H')
  
  ## plotting cell cycle phase ####
  
  p_[[4]] = DimPlot(s_, group.by = 'Phase')+
            theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 15, face = 'bold'))+
            labs(title = 'Cell cycle phase', color = 'Phase')
  
  ## plotting cell plurality ####
  
  p_[[5]] = DimPlot(s_, group.by = 'doublet', order = c('Doublet','Singlet')) +
            theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 15, face = 'bold'))+
            labs(color = 'Plurality', title = 'Doublets (before)')
  
  p_[[6]] = DimPlot(s_, group.by = 'doublet_corrected', order = c('Doublet','Singlet')) +
            theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 15, face = 'bold'))+
            labs(color = 'Plurality', title = 'Doublets (after)')
  
  ## Plotting malignant cells ####
  
  p_[[7]] = DimPlot(s_, cells.highlight = cells_, sizes.highlight = .4)+
            theme(plot.title = element_text(hjust = .5, face = 'bold'), plot.subtitle = element_text(hjust = .5, face = 'bold'),
                  axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
            labs(title = 'Malignant cells')
  
  pdf(file = paste0(sample_,'_preprocessing.pdf'), width = 8, height = 7.5)
    lapply(X = p_, FUN = function(p_){ plot(p_); return(NULL)})
  graphics.off()

  return(NULL)
}

# Preprocessing function ####

# data_ = NULL,             sparse matrix of raw gene counts
# DOUBLETS = NULL,          a named character vector of doublet classification for each cell; if NULL, doublet identification will be performed
# min_features = 200,       minimum number of genes detected in a cell
# max_features = Inf,       maximum number of genes detected in a cell
# max_percent_mt = 10)      maximum percent of mitochondrial genes in a cell
preprocess = function(mat_ = NULL, features_ = NULL, barcodes_ = NULL, sample_ = 'sample', DOUBLETS = NULL, min_features = 200, max_features = Inf, max_percent_mt = 10)      # maximum percent of mitochondrial genes in a cell
{
  message('\n<< Preprocessing sample ', sample_,' will take a while ... >>')
  
  data_ = readMM(file = mat_)                 # reading in count matrix
  features_ = readLines(con = features_)      # reading in gene names
  rownames(data_) = features_
  barcodes_ = readLines(con = barcodes_)      # reading in cell barcodes
  colnames(data_) = barcodes_
  
  message('dimensions of count matrix:',
          '\n number of cells: ', length(colnames(data_)),
          '\n number of genes: ', length(rownames(data_)),
          '\n length of DOUBLETS: ', length(DOUBLETS),
          '\n min features: ', min_features,
          '\n max features: ', max_features,
          '\n max percent MT: ', max_percent_mt)

  s_obj = CreateSeuratObject(counts = data_, min.cells = 10, min.features = 0, project = 'sample')      # creating Seurat object
  
  # QC ####
  
  s_obj[["percent.mt"]] = PercentageFeatureSet(s_obj, pattern = "^MT-")
  s_obj[["percent.mt"]][ is.na(s_obj[["percent.mt"]]) ] = 0
  s_obj = subset(s_obj, subset = min_features <= nFeature_RNA & nFeature_RNA <= max_features & percent.mt <= max_percent_mt)
  
  # Preprocessing ####
  message('\nNormaling the data')
  
  s_obj = NormalizeData(s_obj)                                                  # Library normalization ####
  s_obj = FindVariableFeatures(s_obj)                                           # Finding variable genes ####
  s_obj = ScaleData(s_obj)                                                      # Z transformation ####
  s_obj = RunPCA(s_obj, npcs = 30, ndims.print = 1:2, nfeatures.print = 5)      # PCA dimension reduction ####
  s_obj = RunUMAP(s_obj, dims = 1:30)                                           # UMAP projection ####
  
  # Identifying cell cycle phase ####
  message('\nIdentifying cell cycle phase')
  
  s_obj = CellCycleScoring(s_obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T)
  
  # Identifying doublet ####
  message('\nIdentifying doublets')
  
  if(is.null(DOUBLETS))
  {
    DOUBLETS = doubletIndentify(s_obj)
  }else
  {
    DOUBLETS = DOUBLETS[colnames(s_obj)]
  }
  s_obj$doublet = DOUBLETS
  s_obj$doublet_corrected = s_obj$doublet
  s_obj$doublet_corrected[s_obj$doublet %in% 'Doublet' & ! s_obj$Phase %in% 'G1'] = 'Singlet'     # G1 is senescence/quiescence phase; doublets in proliferative phases (S and G2M) are false positives
  s_obj = subset(s_obj, subset = doublet_corrected %in% 'Singlet')                                # removing doublets
  
  # Identifying malignant cells ####
  message('Identifying malignant cells')
  
  malignant_genes = readLines(con = 'data/malignancy_signature.txt')        # malignant signature
  malignant_cells = malignant_cells_identifier(s_obj, malignant_genes)      # computing p-values for each cell
  
  # f (this part has been made hidden until publication) ####
  # returns a feature vector for current example
  f_vect = f(...)
  
  # Plotting ####
  
  null_ = plotting(s_obj, malignant_cells, sample_)
    
  return(t(f_vect))
}

