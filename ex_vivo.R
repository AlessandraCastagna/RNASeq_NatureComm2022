### Differential expression analyses 
### Analysis in PAIRED

### EdgeR
### DESeq2

# Directory
    #setwd("")

# librerie

    library("DESeq2")
    library("limma")
    library("edgeR")
    library("pheatmap")
    library("gplots")
    library("AnnotationDbi")
    library("org.Hs.eg.db")
    library("qusage")
    library("EGSEA")


## Output folders

                            folder_edgeR_txt = "2019_05_06_Confronti EdgeR/txt/"   #Analisi EdgeR

                            folder_DESeq2_txt = "2019_05_06_Confronti DESeq2/txt/"   #Analisi DESeq2

                            folder_heatmap_edger = "2019_05_06_Heatmap EdgeR/"   #Analisi DESeq2
                            
                            folder_arricchimenti_txt = "2019_05_06_Arricchimenti di espressione genica/txt/"

##            Importo i file
    table <- read.table("raw_counts_ex_vivo.txt", sep="\t", header=T, row.names=1)
    sample.annotation<-read.table("sample_file_ES_AC.txt",header=T, sep="\t")
    rownames(sample.annotation) <- sample.annotation$Label

# Tolgo il replicato problematico

    table = table[, -match( "F2_4N_1", colnames(table))]
    
    sample.annotation = sample.annotation[-match( "F2_4N_1", rownames(sample.annotation)),]
    
## Hallmarks
    H = read.gmt("h.all.v6.2.entrez.gmt")
    gs_annot = data.frame(ID=names(H), GeneSet=names(H))

# C7
    C7 = read.gmt("c7.all.v6.2.entrez.gmt")
    gs_annot_C7 = data.frame(ID=names(C7), GeneSet=names(C7))

# BIOCARTA
    BIOCARTA = read.gmt("c2.cp.biocarta.v6.2.entrez.gmt")
    gs_annot_BIOCARTA = data.frame(ID=names(BIOCARTA), GeneSet=names(BIOCARTA))

# REACTOME
    REACTOME = read.gmt("c2.cp.reactome.v6.2.entrez.gmt")
    gs_annot_REACTOME = data.frame(ID=names(REACTOME), GeneSet=names(REACTOME))

## Design in PAIRED

    Group <- sample.annotation$Group3
    Replicate<- sample.annotation$Replicate
    design <- model.matrix(~0+Group+Replicate)
    colnames(design) <- c(levels(Group), "Replicate")
    
    datamatrix <- table[,c(2:(dim(table)[2]))]
    gene_info <- table[,1]
    gene <- data.frame("ENSEMBL"=row.names(table), "SYMBOL"=table[,1])
    row.names(gene) = gene$ENSEMBL

## Design in PAIRED
    
    Group <- sample.annotation$Colture
    Replicate<- sample.annotation$Replicate
    design <- model.matrix(~0+Group+Replicate)
    colnames(design) <- c(levels(Group), "Replicate")
    
    datamatrix <- table[,c(2:(dim(table)[2]))]
    gene_info <- table[,1]
    gene <- data.frame("ENSEMBL"=row.names(table), "SYMBOL"=table[,1])
    row.names(gene) = gene$ENSEMBL
    
    
    
# ENTREZ ID 
    entrez <- mapIds(org.Hs.eg.db,
                     keys = row.names(gene),
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
    gene$ENTREZ = entrez

# Definisco i contrasti
## Validi sia per EdgeR che per DESeq2

    
    ## confronto cocoltura e singola coltura
    
    
my.contrasts <- makeContrasts(
  
        # M0
        M0_cc_24N - M0_singole_24N,
        M0_cc_4H - M0_singole_4H,
        M0_cc_24H - M0_singole_24H,
        
        # F0
        F0_cc_24N - F0_singole_24N,
        F0_cc_4H - F0_singole_4H,
        F0_cc_24H - F0_singole_24H,
  
        #M1
        M1_cc_4N - M0_cc_24N,
        M1_singole_4N - M0_singole_24N,
        M1_cc_4N - M1_singole_4N,
        ( M1_cc_4N - M0_cc_24N) - (M1_singole_4N - M0_singole_24N),
        
        M1_cc_24N - M0_cc_24N,
        M1_singole_24N - M0_singole_24N,
        M1_cc_24N - M1_singole_24N,
        ( M1_cc_24N - M0_cc_24N) - (M1_singole_24N - M0_singole_24N),
        
        M1_cc_4H - M0_cc_4H,
        M1_singole_4H - M0_singole_4H,
        M1_cc_4H - M1_singole_4H,
        ( M1_cc_4H - M0_cc_4H) - (M1_singole_4H - M0_singole_4H),
        
        M1_cc_24H - M0_cc_24H,
        M1_singole_24H - M0_singole_24H,
        M1_cc_24H - M1_singole_24H,
        ( M1_cc_24H - M0_cc_24H) - (M1_singole_24H - M0_singole_24H),
        
        #M2
        M2_cc_4N - M0_cc_24N,
        M2_singole_4N - M0_singole_24N,
        M2_cc_4N - M2_singole_4N,
        ( M2_cc_4N - M0_cc_24N) - (M2_singole_4N - M0_singole_24N),
        
        M2_cc_24N - M0_cc_24N,
        M2_singole_24N - M0_singole_24N,
        M2_cc_24N - M2_singole_24N,
        ( M2_cc_24N - M0_cc_24N) - (M2_singole_24N - M0_singole_24N),
        
        M2_cc_4H - M0_cc_4H,
        M2_singole_4H - M0_singole_4H,
        M2_cc_4H - M2_singole_4H,
        ( M2_cc_4H - M0_cc_4H) - (M2_singole_4H - M0_singole_4H),
        
        M2_cc_24H - M0_cc_24H,
        M2_singole_24H - M0_singole_24H,
        M2_cc_24H - M2_singole_24H,
        ( M2_cc_24H - M0_cc_24H) - (M2_singole_24H - M0_singole_24H),
        
        #F1
        F1_cc_4N - F0_cc_24N,
        F1_singole_4N - F0_singole_24N,
        F1_cc_4N - F1_singole_4N,
        ( F1_cc_4N - F0_cc_24N) - (F1_singole_4N - F0_singole_24N),
        
        F1_cc_24N - F0_cc_24N,
        F1_singole_24N - F0_singole_24N,
        F1_cc_24N - F1_singole_24N,
        ( F1_cc_24N - F0_cc_24N) - (F1_singole_24N - F0_singole_24N),
        
        F1_cc_4H - F0_cc_4H,
        F1_singole_4H - F0_singole_4H,
        F1_cc_4H - F1_singole_4H,
        ( F1_cc_4H - F0_cc_4H) - (F1_singole_4H - F0_singole_4H),
        
        F1_cc_24H - F0_cc_24H,
        F1_singole_24H - F0_singole_24H,
        F1_cc_24H - F1_singole_24H,
        ( F1_cc_24H - F0_cc_24H) - (F1_singole_24H - F0_singole_24H),
        
        #F2
        F2_cc_4N - F0_cc_24N,
        F2_singole_4N - F0_singole_24N,
        F2_cc_4N - F2_singole_4N,
        ( F2_cc_4N - F0_cc_24N) - (F2_singole_4N - F0_singole_24N),
        
        F2_cc_24N - F0_cc_24N,
        F2_singole_24N - F0_singole_24N,
        F2_cc_24N - F2_singole_24N,
        ( F2_cc_24N - F0_cc_24N) - (F2_singole_24N - F0_singole_24N),
        
        F2_cc_4H - F0_cc_4H,
        F2_singole_4H - F0_singole_4H,
        F2_cc_4H - F2_singole_4H,
        ( F2_cc_4H - F0_cc_4H) - (F2_singole_4H - F0_singole_4H),
        
        F2_cc_24H - F0_cc_24H,
        F2_singole_24H - F0_singole_24H,
        F2_cc_24H - F2_singole_24H,
        ( F2_cc_24H - F0_cc_24H) - (F2_singole_24H - F0_singole_24H),
        
        
        levels=make.names(colnames(design))
)







## Confronto ipossia e normossia 

## confronto cocoltura e singola coltura


my.contrasts <- makeContrasts(
  
        # M0
        M0_cc_4H - M0_cc_24N,
        M0_singole_4H - M0_singole_24N,
        M0_cc_24H - M0_cc_24N,
        M0_singole_24H - M0_singole_24N,
        
        # M1
        M1_cc_4H - M1_cc_4N,
        M1_singole_4H - M1_singole_4N,
        M1_cc_24H - M1_cc_24N,
        M1_singole_24H - M1_singole_24N,
        
        # M2
        M2_cc_4H - M2_cc_4N,
        M2_singole_4H - M2_singole_4N,
        M2_cc_24H - M2_cc_24N,
        M2_singole_24H - M2_singole_24N,
        
        # F0
        F0_cc_4H - F0_cc_24N,
        F0_singole_4H - F0_singole_24N,
        F0_cc_24H - F0_cc_24N,
        F0_singole_24H - F0_singole_24N,
        
        # F1
        F1_cc_4H - F1_cc_4N,
        F1_singole_4H - F1_singole_4N,
        F1_cc_24H - F1_cc_24N,
        F1_singole_24H - F1_singole_24N,
        
        # F2
        F2_cc_4H - F2_cc_4N,
        F2_singole_4H - F2_singole_4N,
        F2_cc_24H - F2_cc_24N,
        F2_singole_24H - F2_singole_24N,
        
        
        levels=make.names(colnames(design))
)




my.contrasts <- makeContrasts(
  
  M1_singole_4H - M1_singole_4N,
  M2_singole_4H - M2_singole_4N,
  F1_singole_4H - F1_singole_4N,
  F2_singole_4H - F2_singole_4N,
  
  levels=make.names(colnames(design))
)


### confronto fibroblasti vs macrofagi


my.contrasts <- makeContrasts(
  
      F0_singole_24N - M0_singole_24N,
      F0_cc_24N - M0_cc_24N,
      F0_singole_4H - M0_singole_4H,
      F0_cc_4H - M0_cc_4H,
      F0_singole_24H - M0_singole_24H,
      F0_cc_24H - M0_cc_24H,
      
      F1_singole_4N - M1_singole_4N,
      F1_cc_4N - M1_cc_4N,
      F1_singole_24N - M1_singole_24N,
      F1_cc_24N - M1_cc_24N,
      F1_singole_4H - M1_singole_4H,
      F1_cc_4H - M1_cc_4H,
      F1_singole_24H - M1_singole_24H,
      F1_cc_24H - M1_cc_24H,
      
      F2_singole_4N - M2_singole_4N,
      F2_cc_4N - M2_cc_4N,
      F2_singole_24N - M2_singole_24N,
      F2_cc_24N - M2_cc_24N,
      F2_singole_4H - M2_singole_4H,
      F2_cc_4H - M2_cc_4H,
      F2_singole_24H - M2_singole_24H,
      F2_cc_24H - M2_cc_24H,

  
  levels=make.names(colnames(design))
)


### confronto singole vs cc


my.contrasts <- makeContrasts(
  
  singole - cc,
  
  levels=make.names(colnames(design))
)


######################
## Analisi con EdgeR##
######################


##        Creo oggetto DGEList
    y <- DGEList(counts=datamatrix,genes=gene_info)

##   Filtering
    keep <- rowSums(cpm(y)>1) >= 3      #un gene deve apparire in almeno 3 campioni
    y <- y[keep, , keep.lib.sizes=FALSE]         #ricalcola la dimensione della matrice
    #summary(cpm(y))

##             Normalization
    y <- calcNormFactors(y)
    y<-estimateCommonDisp(y)
    y<-estimateTagwiseDisp(y)
##                   Stima della dispersione
    y <- estimateDisp(y, design)
    
    fit <- glmQLFit(y, design)

## Estimating the dispersion
    plotBCV(y)

    plotQLDisp(fit)
############################


#### Ciclo per ogni contrasto definito
## Viene calcolata l'espressione differenziale
## Correlazione
## Gene Ontology
## Kegg pathways

for (xxx in 1:(dim(my.contrasts)[2])) {
            print("Confronto EDGER n°")  
            print(xxx)
            
            print(colnames(my.contrasts)[xxx])
            
            lista = my.contrasts [ which (my.contrasts [,xxx] != 0 ) , xxx]
            lista = rownames(as.matrix(lista))
            
            
            
            if (length(lista) == 2) { sample.annotation.selected = sample.annotation[which(sample.annotation$Group3 == lista[1] 
                                                                                           | sample.annotation$Group3 == lista[2] ), ] 
                                     
            }
            
            if (length(lista) == 4) { sample.annotation.selected = sample.annotation[which(sample.annotation$Group3 == lista[1] 
                                                                                           | sample.annotation$Group3 == lista[2]
                                                                                           | sample.annotation$Group3 == lista[3]
                                                                                           | sample.annotation$Group3 == lista[4]), ] 
                                      
            }
            
            main = colnames(my.contrasts)[xxx]
            main = gsub("-", "VS", main, "_EdgeR_")
            
            sample.annotation.selected = sample.annotation.selected[order(sample.annotation.selected$Replicate),]
            
            # Controllo com'è composta la matrice di annotazioni
            # Avendo tolto un replicato, bisogna controllare che la selezione sia fatta in modo opportuno
            # Infatti si tratta di un'analisi in PAIRED
            
            if ( dim(sample.annotation.selected)[1] %% 2 == 1){
              
              sample.annotation.selected = sample.annotation.selected[ which(sample.annotation.selected$Replicate == 2 | sample.annotation.selected$Replicate == 3 ) ,]
              
            }
            
            # Analisi differenziale
            
            print("Analisi differenziale")
            
            qlf <- glmQLFTest(fit, contrast=my.contrasts[,xxx])
            res<-topTags(qlf,n=30000)
            
            nome.file = paste0(main, "_PAIRED_")
            
            write.table(res, paste0(folder_edgeR_txt, nome.file,"_differential_expression.txt"), sep="\t",quote = F, col.names = NA)
            write.table(res, paste0(folder_edgeR_punto, nome.file,"_differential_expression.xls"), sep="\t",quote = F, col.names = NA)
            write.table(res, paste0(folder_edgeR_virgola, nome.file,"_differential_expression.xls"), dec = ",", sep="\t",quote = F, col.names = NA)
            
            # Geni differenzialmente espressi significativi
            gene.selected <- rownames(res[res$table$FDR<0.05,])
            
            print("Geni differenzialmente espressi significativi con FDR < 0,05")
            print(length(gene.selected))
            
            # Correlazione campioni
           
            print("Correlazione campioni")
            
            data <- y$pseudo.counts
            
            mean = c()
            
            for (pluto in 1:dim(data)[1]){
              media = mean(data[pluto,])
              mean = c(mean, media)
            }
            
            valore = 50
            
            select=data[which((mean)>=valore),rownames(sample.annotation.selected)]
            
            matrice <- matrix(nrow=dim(select)[2], ncol=dim(select)[2], data = NA)
            colnames(matrice) <- colnames(select)
            rownames(matrice) <- colnames(select)
            
            for (i in 1:(dim(select)[2]-1)){
              matrice[i,i]=1
              for(j in (i+1):(dim(select)[2])){
                alpha <- cor.test(select[,i], select[,j], method = "pearson")
                matrice[i,j]=alpha$estimate
                matrice[j,i]=alpha$estimate
              }
            }
            
            matrice[dim(matrice)[1],dim(matrice)[2]]=1
            
            write.table(matrice, paste0(folder_heatmap_edger, nome.file,"_Correlazione.txt"), sep="\t", col.names = NA, quote =F)
            
            print("Heatmap Correlazione")
            x = pheatmap(matrice, 
                         scale = "none", 
                         fontsize = 8,
                         legend = TRUE, 
                         border_color = NA,
                         cellwidth = 20, cellheight = 20,
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         main = "Correlazione",
                         filename = paste0(folder_heatmap_edger, nome.file,"_Corr_heatmap.png"))
            
            
            if ( length ( gene.selected ) > 0 ) {
                  
                  print("C'è almeno un gene significativo")
                    
                  # tabella di espressione
                  
                  table.selected <- y$pseudo.counts[gene.selected,rownames(sample.annotation.selected)]
                  rownames(table.selected) <- gene[rownames(table.selected),"SYMBOL"]
                  
                  # heatmap livelli di espressione dei campioni coinvolti nel confronto
                  
                  print("Salva Heatmap livelli di espressione")
                  
                  x = pheatmap(table.selected, 
                               scale = "row", 
                               show_rownames = FALSE,
                               border_color = NA,
                               cluster_cols = FALSE,
                               cluster_rows = TRUE,
                               clustering_distance_rows = "euclidean",
                               main = main,
                               filename = paste0(folder_heatmap_edger, nome.file,"_heatmap_espressione.png"))
                  
                  #geni significativi positivi
                  gene.selected.2 <- rownames(res[which ( (res$table$logFC) > 0),])
                  gene.selected.positive = intersect(gene.selected,gene.selected.2)
                  print("Geni significativi positivi")
                  print(length(gene.selected.positive))
                  
                  #geni significativi negativi
                  gene.selected.3 <- rownames(res[which ( (res$table$logFC) < 0),])
                  gene.selected.negative = intersect(gene.selected,gene.selected.3)
                  print("Geni significativi negativi")
                  print(length(gene.selected.negative))
                  
                  #geni con |logFC|>1
                  gene.selected.4 <- rownames(res[which ( (abs(res$table$logFC) > 1)),])
                  gene.selected.logFC = intersect(gene.selected,gene.selected.4)
                  print("Geni con |logFC| > 1")
                  print(length(gene.selected.logFC))
                  
                 
                  
                     
                  
            }
            
              
              
              # tabella di espressione
              print("Salva tabella di espressione")
              table.selected <- y$pseudo.counts[,rownames(sample.annotation.selected)]
              rownames(table.selected) <- gene[rownames(table.selected),"SYMBOL"]
              
              write.table(table.selected, paste0(folder_edgeR_txt,nome.file,"_espressione.txt"), sep="\t",quote = F, col.names = NA)
              }
    
    
