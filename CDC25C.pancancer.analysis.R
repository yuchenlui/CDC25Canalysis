#!/usr/bin/env Rscript
# =====================================================================
# CDC25C.pancancer.analysis.R - CDC25C single-gene pan-cancer pipeline (single file)
# 2026-10-12
#
# [Data] Expression input, single source (~9.5 GB):
#   UCSC Xena Toil Recomputed pan-cancer TPM matrix
#   Download page: https://xenabrowser.net/datapages/?cohort=TCGA%20Target-GTEx
#   dataset "TcgaTargetGtex_rsem_gene_tpm" (Hugo symbol version,
#   file name TcgaTargetGtex_rsem_gene_tpm.symbol.tsv)
#   place at data/TcgaTargetGtex_rsem_gene_tpm.symbol.tsv under the project root
#
# [Other inputs] per-step reference data (bundled directories):
#   07.preSurvival/ (Xena TCGA survival data, steps 08-12),13.cliCor/ 15.TMBcor/ 17.MSIcor/ 20.estimateCor/
#   21.CIBERSORT/ 23.geneCor/ 25.GSEA/ database/, 29.../symbol.txt (expression matrix for drug sensitivity),
#   33.IMvigor210/IMvigor210/ (immunotherapy cohort)
#
# [Output]
#   data/CDC25C.tpm.denseDataOnlyDownload.tsv (input to 06.diff)
#   data/CDC25CexprTime.txt (expression input to steps 08-29, TCGA samples only)
#   <step dir>/<output dir>/ per-step results; report/CDC25C/ collected copy
#
# [Usage] from the project root: Rscript CDC25C.pancancer.analysis.R
# [Notes] Expression scale: log2(TPM+1) (Toil); High/Low split at the median,
#         equivalent to the FPKM-based grouping; each step runs in tryCatch,
#         a failed step is logged and does not stop the pipeline
# =====================================================================

suppressPackageStartupMessages({
  .pkgs <- c("data.table", "ggplot2", "dplyr", "pheatmap", "RColorBrewer", "survival",
             "survminer", "stringr", "xlsx", "forestplot", "ggpubr", "fmsb", "patchwork",
             "reshape2", "tidyr", "scales", "openxlsx")
  .miss <- .pkgs[!vapply(.pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss)) {
    message("installing missing CRAN packages: ", paste(.miss, collapse = ", "))
    try(install.packages(.miss, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"), silent = TRUE)
  }
  .bioc <- c("org.Hs.eg.db", "clusterProfiler", "limma", "edgeR", "GSVA", "ggplotify",
             "GSEABase", "fgsea", "enrichplot", "DOSE", "preprocessCore", "survcomp")
  .miss2 <- .bioc[!vapply(.bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss2)) {
    message("installing missing Bioconductor packages: ", paste(.miss2, collapse = ", "))
    try({
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(.miss2, ask = FALSE, update = FALSE)
    }, silent = TRUE)
  }
})

# ------------------------------ global parameters ------------------------------
gene <- "CDC25C"
ParentWD <- normalizePath(".")
cohort_focus <- c("PAAD", "COAD", "LIHC", "STAD")
assign("gene", gene, envir = globalenv())
assign("cohort_to_analysis", cohort_focus, envir = globalenv())
assign("significant_cohort", c(), envir = globalenv())


# =====================================================================
# 1. Expression extraction (single source: Toil matrix)
# =====================================================================
dataDir <- file.path(ParentWD, "data")
toil_big <- file.path(dataDir, "TcgaTargetGtex_rsem_gene_tpm.symbol.tsv")
stopifnot(file.exists(toil_big))

scan_genes <- function(path, want, chunk = 400L) {
  con <- file(path, "r"); on.exit(close(con))
  hdr <- readLines(con, n = 1L)
  samples <- strsplit(hdr, "\t", fixed = TRUE)[[1]][-1L]
  ns <- length(samples)
  mat <- matrix(NA_real_, nrow = ns, ncol = length(want), dimnames = list(NULL, want))
  found <- rep(FALSE, length(want)); nb <- 0L
  repeat {
    ln <- readLines(con, n = chunk)
    if (length(ln) == 0L) break
    nb <- nb + length(ln)
    gn <- sub("\t.*$", "", substr(ln, 1L, 80L))
    for (k in seq_along(want)) {
      if (found[k]) next
      i <- match(want[k], gn)
      if (!is.na(i)) {
        v <- strsplit(ln[i], "\t", fixed = TRUE)[[1]][-1L]
        if (length(v) != ns) stop("column count mismatch: ", basename(path))
        mat[, k] <- as.numeric(v); found[k] <- TRUE
      }
    }
    if (nb %% 20000L == 0L) message("    scanned ", nb, " lines ...")
    if (all(found)) break
  }
  if (any(!found)) stop("gene(s) not found in the matrix: ", paste(want[!found], collapse = ","))
  list(samples = samples, mat = mat)
}

message("[extract] scanning Toil matrix (~10 GB, takes a few minutes) ...")
toil <- scan_genes(toil_big, gene)
sid <- toil$samples
val <- toil$mat[, gene]        # log2(TPM+1)

## ---- Part A: Xena denseDataOnlyDownload-equivalent format (input to 06.diff) ----
val_xena <- round(log2(pmax(2^val - 0.99, 1e-9)), 3)   # to log2(TPM+0.01) scale
study <- ifelse(grepl("^TCGA", sid), "TCGA", ifelse(grepl("^GTEX", sid), "GTEX", "OTHER"))
tc <- suppressWarnings(as.integer(substr(sid, 14L, 15L)))
stype <- ifelse(study == "GTEX", "Normal Tissue",
         ifelse(is.na(tc), "Primary Tumor",
         ifelse(tc == 2L, "Recurrent Tumor",
         ifelse(tc %in% c(6L, 7L), "Metastatic",
         ifelse(tc >= 10L & tc <= 19L, "Solid Tissue Normal",
         ifelse(tc >= 20L, "Control Analyte", "Primary Tumor"))))))
fwrite(data.table(sample = sid, samples = sid, setNames(list(val_xena), gene),
                  `_sample_type` = stype, `_study` = study),
       file.path(dataDir, paste0(gene, ".tpm.denseDataOnlyDownload.tsv")),
       sep = "\t", quote = FALSE)
message("[extract] Part A done -> data/", gene, ".tpm.denseDataOnlyDownload.tsv")

## ---- Part B: exprTime (TCGA samples only, input to steps 08-29) ----
is_tcga <- grepl("^TCGA", sid)
tsid <- sid[is_tcga]; tval <- val[is_tcga]
cancer <- vapply(strsplit(tsid, "-"), function(x) if (length(x) >= 2L) x[2] else NA_character_, character(1))
code <- vapply(strsplit(tsid, "-"), function(x) if (length(x) >= 4L) x[4] else NA_character_, character(1))
ttype <- ifelse(substr(code, 1L, 1L) == "1", "Normal", "Tumor")
ok <- !is.na(cancer) & nzchar(cancer)
fwrite(data.table(ID = tsid[ok], Type = ttype[ok], cancer_type = cancer[ok],
                  setNames(list(tval[ok]), gene)),
       file.path(dataDir, paste0(gene, "exprTime.txt")), sep = ",", quote = FALSE)
message("[extract] Part B done -> data/", gene, "exprTime.txt (", sum(ok), " TCGA samples, ",
        length(unique(cancer[ok])), " cancer types)")

# =====================================================================
# 2. Report collection
# =====================================================================
cp_report <- function(dir) {
  from <- file.path(ParentWD, dir)
  to <- file.path(ParentWD, "report", gene)
  if (!dir.exists(to)) dir.create(to, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(from)) { warning("[report] directory not found, skipped: ", from); return(invisible(FALSE)) }
  invisible(file.copy(from, to, recursive = TRUE, copy.date = TRUE))
}

# Step runner: setwd to <step dir>/<output dir>, then eval the inlined step code
run_step <- function(step_dir, outdir, codes) {
  odir <- file.path(ParentWD, step_dir, outdir)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  owd <- setwd(odir)
  on.exit(setwd(owd), add = TRUE)
  invisible(TRUE)
}

# =====================================================================
# 3. Analysis steps (inlined from each step script; executed inside <step dir>/<output dir>)
# =====================================================================

step_06_diff <- function() {
  outdir <- paste0("1.", gene, ".diff.boxplot")
  message("[step_06_diff] start")
  setwd(file.path(ParentWD, "06.diff"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(dplyr)


  dat=read.table(paste0(ParentWD,"/data/",gene,".tpm.denseDataOnlyDownload.tsv"),sep="\t",header = T)
  pandata=read.csv(paste0(ParentWD,"/06.diff/pancancer.csv"))
  cancerlist=unique(pandata$cancer_type)
  for(i in cancerlist){
    cancer=pandata%>%filter(cancer_type==i)
    p=merge(cancer,dat,by="samples")
    write.csv(p,paste0("csv/",i,".csv"))
  }

  filename=list.files("csv")
  filename=paste0("csv/",filename)
  merge.df=read.csv(filename[1],header = T)
  for (i in 2:length(filename)){
    new.df=read.csv(filename[i],header = T)
    merge.df=rbind(merge.df,new.df)
    }
  merge.df$sample_type=ifelse(merge.df$X_study=="TCGA","Tumor","Normal")
  write.csv(merge.df,paste0(gene,".merge.csv"),row.names = F)

  dat=merge.df
  colnames(dat)[5]=gene
  dat[,5]=round(log2(2^dat[,5]-0.01+1),3)
  dat=dat[which(dat$X_sample_type!="Solid Tissue Normal"),]
  library(ggpubr)
  p=ggboxplot(dat, # 'dat' from the environment as plot data
              x="cancer_type", y=gene, color = "sample_type", 
              ylab=paste0(gene," expression"),xlab="",
              palette = c("#244B9E", "firebrick"),ggtheme = theme_bw())+rotate_x_text(60)+stat_compare_means(aes(group=sample_type),
                       method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                       label = "p.signif")+theme(legend.title = element_blank(),axis.text.x = element_text(size = 12,colour = "black"),axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12,colour = "black"),legend.text = element_text(size = 12,colour = "black"))

  pdf(file=paste0(gene,"_boxplot_tpm_UCSC_TCGA_GTEx.pdf"),width=12,height=6)  
  print(p)
  dev.off()

  unlink("csv", recursive=TRUE)


  if(!dir.exists("stage")){ dir.create(paste0("stage"))
    setwd("stage")
  }else{message("There are some reports, please check out!")
    setwd("stage")}

  rt=dat#genename.merge.csv
  rt$cancer_type=factor(rt$cancer_type)
  gtex=rt[grep("^GTEX",rt$samples),]
  colnames(gtex)[5]="gene"
  rt1=rt[grep("^TCGA",rt$samples),]

  rownames(rt1)=rt1$sample
  cli=read.table(paste0(ParentWD,"/13.cliCor/clinical.txt"),sep="\t",header=T,check.names=F,row.names=1) 
  cli$samples=rownames(cli)
  clinical=colnames(cli)[1]


  if(nchar(row.names(rt1)[1])!=nchar(row.names(cli)[1])){
    row.names(rt1)=gsub(".$","",row.names(data))}

  sameSample=intersect(row.names(rt1),row.names(cli))
  sameData=rt1[sameSample,]
  sameClinical=cli[sameSample,]
  cliExpData=merge(sameClinical,sameData,by="samples")
  colnames(cliExpData)[2]="sameClinical"
  colnames(cliExpData)[6]="gene"
  rt1=cliExpData%>%dplyr::select(c("sameClinical","samples","gene","cancer_type"))
  rt1$gene=as.numeric(as.character(rt1$gene))
  gtex1=data.frame(sameClinical=rep("Normal",nrow(gtex)),samples=as.character(gtex$samples),gene=as.numeric(gtex$gene),cancer_type=as.character(gtex$cancer_type))
  rt1$sameClinical=as.character(rt1$sameClinical)
  rt1$gene=as.numeric(as.character(rt1$gene))
  rt1$samples=as.character(rt1$samples)
  rt1$cancer_type=as.character(rt1$cancer_type)
  library(ggpubr)
  library(dplyr)
  cliExpData=rbind(rt1,gtex1)
  cliExpData$sameClinical=factor(cliExpData$sameClinical,levels =c("Normal","Stage I","Stage II","Stage III","Stage IV"),labels = c("Normal","Stage I","Stage II","Stage III","Stage IV"))

  for(i in levels(rt[,"cancer_type"])){
    cliExpData1=cliExpData[(cliExpData[,"cancer_type"]==i),]
    if(nrow(cliExpData1)==0){next}
    cliExpData1$sameClinical=factor(cliExpData1$sameClinical)
    group=levels(cliExpData1$sameClinical)
    if(length(group)==1){next}
    comp=combn(group,2)
    my_comparisons=list()
    for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

    boxplot=ggboxplot(cliExpData1, x="sameClinical", y="gene", color="sameClinical",
                      xlab=clinical,
                      ylab=paste(gene,"expression"),
                      legend.title=clinical,
                      title=paste0("Cancer: ",i),
                      add = "jitter")+ 
      stat_compare_means(comparisons = my_comparisons)
    pdf(file=paste0(gene,".",clinical,".Normal.",i,".pdf"),width=5.5,height=5)
    print(boxplot)
    dev.off()
  }
  setwd("..")


  library(ggpubr)
  library(dplyr)
  if(!dir.exists("paired")){ dir.create(paste0("paired"))
    setwd("paired")
  }else{message("There are some reports, please check out!")
    setwd("paired")}

  merge.df=read.csv(paste0("../",gene,".merge.csv"))

  colnames(merge.df)[5]="gene"
  df=as.data.frame(merge.df)%>%dplyr::filter(X_study=="TCGA")%>%dplyr::select(c("sample","cancer_type","gene","X_sample_type"))
  df$sample=as.character(df$sample)
  df$person=stringr::str_sub(df$sample,1,-4)
  df$X_sample_type=as.character(df$X_sample_type)
  df$cancer_type=factor(df$cancer_type)
  for( i in levels(df$cancer_type)){

  df1=df[(df[,"cancer_type"]==i),]
  df_normal=df1%>%dplyr::filter(X_sample_type=="Solid Tissue Normal")
  if(nrow(df_normal)==0){
    print("no normal")
    next
    }else{print("ok")}



  to_del=c()
  df_paired=data.frame()
  for(j in 1:nrow(df_normal)){
    temp1=df[which(df$person==df_normal$person[j]&df$X_sample_type=="Primary Tumor"),]
    if(nrow(temp1)==1){
    df_paired=rbind(df_paired,temp1)
    }else if(nrow(temp1)==0){
      temp2=df[which(df$person==df_normal$person[j]&df$X_sample_type=="Metastatic"),]
      if(nrow(temp2)==0){
        to_del=c(to_del,j)
        next
        }
      df_paired=rbind(df_paired,temp2)
    }
  }
  if(length(to_del)!=0){
  df_normal=df_normal[-to_del,]
  }

  if(nrow(df_normal)<2|nrow(df_paired)<2){
    print(paste("too less samples",nrow(df_normal),nrow(df_paired)))
    next
    }
  df_paired=rbind(df_paired,df_normal)
  if(nrow(df_paired)==0){next}
  df_paired$X_sample_type=ifelse(grepl("Normal",df_paired$X_sample_type),"Normal","Tumor")
  colnames(df_paired)[4]="Type"
  df_paired$Type=factor(df_paired$Type,levels = c("Normal","Tumor"))
  p <- ggpaired(df_paired, x = "Type", y = "gene",
                  fill =  "Type",
                  line.color = "black", line.size = 0.4,
                  ylab = paste0(gene," expression level [log(TPM+1)]"))+
    scale_fill_manual(values = c(
      scales::alpha("#244B9E", 0.8), # alpha 0.8
      scales::alpha("firebrick", 0.8)  # alpha 0.8
    ))+
    stat_compare_means(paired = TRUE,method = "t.test",label="p.signif",label.x = 1.5)+theme(legend.position = "none",axis.title.x = element_blank())
  Sys.sleep(1)
  pvalue=compare_means(as.formula(paste0("gene","~Type")),df_paired,method = "t.test",paired = TRUE)$p
  Sys.sleep(1)

  pdf(file=paste0(gene,".",i,".pvalue.",pvalue,"_paird_boxplot_UCSC.pdf"),width=4,height=5)  
  print(p)
  dev.off()

  Sys.sleep(1)
  }

  setwd("..")
  setwd(ParentWD)
  cp_report(file.path("06.diff", outdir))
  message("[step_06_diff] done")
}

step_08_survival <- function() {
  outdir <- paste0("2.", gene, ".survivalplot")
  message("[step_08_survival] start")
  setwd(file.path(ParentWD, "08.survival"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)


  library(limma)

  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  files=dir(paste0(ParentWD,"/07.preSurvival/"))
  files=grep(".survival.tsv$",files,value=T)
  files=paste0(ParentWD,"/07.preSurvival/",files)

  surTime=data.frame()
  for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
  }
  colnames(surTime)=c("futime","fustat")
  surTime=surTime[stringr::str_sub(row.names(surTime),-1,-1)=="A",]
  rownames(exp)=exp$samples
  sameSample=intersect(row.names(surTime),row.names(exp))
  surTime=surTime[sameSample,]
  exp=exp[sameSample,]
  surData=cbind(surTime,exp)
  surData=cbind(id=row.names(surData),surData)
  surData=surData[stringr::str_sub(row.names(surData),-3,-2)=="01",]

  library(survival)
  library(survminer)

  surPlot=list()
  rt=surData
  rt$futime=rt$futime/365
  gene=colnames(rt)[7]
  pFilter=0.05
  rt[,"cancer_type"]=factor(rt[,"cancer_type"])

  for(i in levels(rt[,"cancer_type"])){
  	rt1=rt[(rt[,"cancer_type"]==i),]
  	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  	pValue=1-pchisq(diff$chisq,df=1)
  	if(pValue<pFilter){
  		if(pValue<0.001){
  			pValue="p<0.001"
  		}else{
  			pValue=paste0("p=",sprintf("%.03f",pValue))
  		}
  		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)

  		surPlot=ggsurvplot(fit, 
  				    data=rt1,
  				    title=paste0("Cancer: ",i),
  				    pval=pValue,
  				    pval.size=6,
  				    legend.labs=c("high","low"),
  				    legend.title=paste0(gene," levels"),
  				    font.legend=12,
  				    xlab="Time(years)",
  				    ylab="Overall survival",
  				    break.time.by = 1,
  				    palette=c("firebrick","#244B9E"),
  				    conf.int=F,
  				    fontsize=4,
  				    risk.table=FALSE)

  		pdf(file=paste0(gene,".survival.OS.",i,".pdf"),onefile = FALSE,
  			    width = 6, 
  			    height =6)             
  		print(surPlot)
  		dev.off()
  	}
  }
  for(i in levels(rt[,"cancer_type"])){
    print(paste0("analysising ",i))
    rt1=rt[(rt[,"cancer_type"]==i),]
  	rt3=rt1
  	rt3$group1=ifelse(rt3[,gene]>quantile(rt3[,gene])[4],"high",ifelse(rt3[,gene]<=quantile(rt3[,gene])[2],"low",NA))
  	rt2=rt3[!is.na(rt3$group1),]
  	diff=survdiff(Surv(futime, fustat) ~group1,data = rt3)
  	pValue=1-pchisq(diff$chisq,df=1)
  	if(pValue<0.08){
  	  if(pValue<0.001){
  	    pValue="p<0.001"
  	  }else{
  	    pValue=paste0("p=",sprintf("%.03f",pValue))
  	  }
  	  fit <- survfit(Surv(futime, fustat) ~ group1, data = rt3)

  	  surPlot=ggsurvplot(fit, 
  	                     data=rt3,
  	                     title=paste0("Cancer: ",i),
  	                     pval=pValue,
  	                     pval.size=6,
  	                     legend.labs=c("high","low"),
  	                     legend.title=paste0(gene," levels"),
  	                     font.legend=12,
  	                     xlab="Time(years)",
  	                     ylab="Overall survival",
  	                     break.time.by = 1,
  	                     palette=c("firebrick","#244B9E"),
  	                     conf.int=F,
  	                     fontsize=4,
  	                     risk.table=FALSE)

  	  pdf(file=paste0(gene,".survival.quantile.OS.",i,".pdf"),onefile = FALSE,
  	      width = 6, 
  	      height =6)             
  	  print(surPlot)
  	  dev.off()
  	}
  }
  setwd(ParentWD)
  cp_report(file.path("08.survival", outdir))
  message("[step_08_survival] done")
}

step_09_cox_os <- function() {
  outdir <- paste0("cox.", gene, ".OS.univariate_cox_forest_plot")
  message("[step_09_cox_os] start")
  setwd(file.path(ParentWD, "09.cox"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)


  library(dplyr)
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  files=dir(paste0(ParentWD,"/07.preSurvival/"))
  files=grep(".survival.tsv$",files,value=T)
  files=paste0(ParentWD,"/07.preSurvival/",files)

  surTime=data.frame()
  for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
  }
  colnames(surTime)=c("futime","fustat")
  surTime=surTime[stringr::str_sub(row.names(surTime),-1,-1)=="A",]
  rownames(exp)=exp$samples
  sameSample=intersect(row.names(surTime),row.names(exp))
  surTime=surTime[sameSample,]
  exp=exp[sameSample,]
  surData=cbind(surTime,exp)
  surData=cbind(id=row.names(surData),surData)
  surData=surData[stringr::str_sub(row.names(surData),-3,-2)=="01",]

  library(survival)
  library(survminer)

  surPlot=list()
  rt=surData
  rt$futime=rt$futime/365
  gene=colnames(rt)[7]
  pFilter=0.05
  rt[,"cancer_type"]=factor(rt[,"cancer_type"])




  library(survival)
  library(forestplot)

  options(scipen = 0)

  outTab=data.frame()
  for(i in levels(rt[,"cancer_type"])){
  	rt1=rt[(rt[,"cancer_type"]==i),]
  	cox=coxph(Surv(futime, fustat) ~rt1[,gene] , data = rt1)
  	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]

    # export for downstream processing
    HR=coxSummary$conf.int[,"exp(coef)"]
    HR.95L=confint(cox)[,'2.5 %']
    HR.95H=confint(cox)[,'97.5 %']

    # clip values to keep points within the axis range
    if(HR > 5){
        HR = 5/(1+exp(-HR))
      if(HR.95L < 0){
        if(HR.95L < logb(0.1)){
          HR.95L = -0.1/(1+exp(HR.95L))
        }
        HR.95H= logb(10)/(1+exp(-HR.95H))
      }else if(HR.95L > 0){
        HR.95L = logb(2)/(1+exp(-HR.95L))
        HR.95H = logb(7)/(1+exp(-HR.95H))
      }else{
        HR.95L = 0
        HR.95H = logb(10)/(1+exp(-HR.95H))
      }
    }else if(HR.95H > 5){
      HR = 5/(1+exp(-HR))
      if(HR.95L < 0){
        if(HR.95L < logb(0.1)){
          HR.95L = -0.1/(1+exp(HR.95L))
        }
        if(HR < 1){
          HR.95L = logb(HR) + HR.95L
        }  
        HR.95H= logb(10)/(1+exp(-HR.95H))
      }else{
        HR.95L = logb(2)/(1+exp(-HR.95L))
        HR.95H = logb(10)/(1+exp(-HR.95H))
      }
    }
    if(HR.95L < logb(0.1)){
      HR.95L = -0.1/(1+exp(HR.95L))
      if(HR < 1){
        HR.95L =  logb(HR)+HR.95L
      }               
      if(HR.95H > 10){
        HR= 5/(1+exp(-HR))
        HR.95H= logb(10)/(1+exp(-HR.95H))
        }

    }
    if(HR.95H > 3){
      HR = 5/(1+exp(-HR))
      HR.95H = logb(10)/(1+exp(-HR.95H))
      if(HR.95L < 0){
        if(HR.95L < logb(0.1)){
          HR.95L = -0.1/(1+exp(HR.95L))
        }
        if(HR < 1){
          HR.95L = logb(HR) + HR.95L
        }  
        HR.95H= logb(10)/(1+exp(-HR.95H))
      }else{
        HR.95L = logb(2)/(1+exp(-HR.95L))
        HR.95H = logb(10)/(1+exp(-HR.95H))
      }
    }
  	outTab=rbind(outTab,
  	             cbind(cancer=i,
  	                   HR=HR,
  	                   HR.95L=exp(HR.95L),
  	                   HR.95H=exp(HR.95H),
  	                   #HR.95L=ifelse(coxSummary$conf.int[,"lower .95"]<0.1,0.1,coxSummary$conf.int[,"lower .95"]),
  	                   #HR.95H=ifelse(coxSummary$conf.int[,"upper .95"]>10,coxSummary$conf.int[,"exp(coef)"]+1,coxSummary$conf.int[,"upper .95"]),
  			           pvalue=format(coxP, scientific=F)))
  }
  outTab$pvalue=as.numeric(as.character(outTab$pvalue))
  outTab1=outTab%>%arrange(pvalue,desc(HR))



  significant_cohort=union(significant_cohort,outTab1[which(outTab1$pvalue<0.05),"cancer"])
  write.table(outTab1,file=paste0(gene,".OS.cox.result.txt"),sep="\t",row.names=F,quote=F) 



  bioForest=function(coxFile=null,forestFile=null,forestCol=null){

  	dff=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  	dff1=dff%>%arrange(pvalue,desc(HR))
  	data=as.data.frame(dff1)
  	HR=data[,1:3]
  	hr=sprintf("%.3f",HR[,"HR"])
  	hrLow=sprintf("%.3f",HR[,"HR.95L"])
  	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  	pVal=data[,"pvalue"]
  	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  	clrs <- fpColors(box=forestCol,line="#e10b2b", summary="#e10b2b")     
  	tabletext <- 
  	  list(c(NA, rownames(HR)),
  	       append("pvalue", pVal),
  	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
  	p=forestplot(tabletext, 
  	             rbind(rep(NA, 3), HR),
  	             col=clrs,
  	             graphwidth=unit(50, "mm"),
  	             xlog=T,
  	             lwd.ci=1,
  	             boxsize=0.6,
  	             xlab="Hazard ratio",
  	             txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))	)

  	pdf(file=forestFile,width = 9,height = 6)
  	print(p)
  	dev.off()
  }


  bioForest(coxFile=paste0(gene,".OS.cox.result.txt"),forestFile=paste0(gene,".OS.forest.pdf"),forestCol="#e10b2b")
  setwd(ParentWD)
  cp_report(file.path("09.cox", outdir))
  message("[step_09_cox_os] done")
}

step_10_cox_dss <- function() {
  outdir <- paste0("DSS.", gene, ".DSS.univariate_cox_forest_plot")
  message("[step_10_cox_dss] start")
  setwd(file.path(ParentWD, "10.DSS"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)
  library(limma)
  library(survival)
  library(survminer)
  library(forestplot)
  # lower/upper limits
  min_limit = 0.001
  max_limit = 10
  # squash extreme values into a bounded range
  transform_value <- function(x, min_limit, max_limit) {
    if(x<max_limit & x>min_limit){return(x)}
    # map to a bounded interval: sigmoid on the log scale
    scaled_value <- 1 / (1 + exp(-log(x)))  # sigmoid
    # map values to [min_limit, max_limit]
    result <- min_limit + (max_limit - min_limit) * scaled_value
    return(result)
  }

  pFilter=0.05                                                      
  cliFile="../Survival_SupplementalTable_S1_20171025_xena_sp"

  if("futime"%in%colnames(rt)){
  rt=rt%>%dplyr::select(-futime,-fustat)
  }
  cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)
  cli=cli[,c("DSS.time","DSS")]
  cli=na.omit(cli)
  colnames(cli)=c("futime","fustat")
  rt$cancer_type=factor(rt$cancer_type)

  outTab=data.frame()
  for(i in levels(rt[,"cancer_type"])){
  	rt1=rt[(rt[,"cancer_type"]==i),]

  	data=cbind(rt1,gene=rt1[,gene])
  	data=as.matrix(data[,c(gene,"gene")])
  	if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
  		row.names(data)=gsub(".$","",row.names(data))}
  	data=avereps(data)
  	sameSample=intersect(row.names(data),row.names(cli))
  	sameData=data[sameSample,]
  	sameCli=cli[sameSample,]
  	rt1=cbind(sameCli,sameData)
  	rt1$futime=rt1$futime/365
  	if(nrow(rt1)<3){next}

  	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  	pValue=1-pchisq(diff$chisq,df=1)
  	if(pValue<pFilter){
  		if(pValue<0.001){
  			pValue="p<0.001"
  		}else{
  			pValue=paste0("p=",sprintf("%.03f",pValue))
  		}
  		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)

  		surPlot=ggsurvplot(fit, 
  				    data=rt1,
  				    title=paste0("Cancer: ",i),
  				    pval=pValue,
  				    pval.size=6,
  				    legend.labs=c("high","low"),
  				    legend.title=paste0(gene," levels"),
  				    font.legend=12,
  				    xlab="Time(years)",
  				    ylab="Disease-specific survival",
  				    break.time.by = 1,
  				    palette=c("red","blue"),
  				    conf.int=F,
  				    fontsize=4,
  				    risk.table=FALSE)
  		pdf(file=paste0(gene,".DSS.",i,".pdf"),onefile = FALSE, width = 6, height =6)
  		print(surPlot)
  		dev.off()
  	}


  	cox=coxph(Surv(futime, fustat) ~ gene, data = rt1)
  	coxSummary = summary(cox)
  	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  	# export for downstream processing
  	HR=coxSummary$conf.int[,"exp(coef)"]
  	HR.95L=confint(cox)[,'2.5 %']
  	HR.95H=confint(cox)[,'97.5 %']

  	# clip values to keep points within the axis range
  	if(HR > 5){
  	  HR = 5/(1+exp(-HR))
  	  if(HR.95L < 0){
  	    if(HR.95L < logb(0.1)){
  	      HR.95L = -0.1/(1+exp(HR.95L))
  	    }
  	    HR.95H= logb(10)/(1+exp(-HR.95H))
  	  }else if(HR.95L > 0){
  	    HR.95L = logb(2)/(1+exp(-HR.95L))
  	    HR.95H = logb(7)/(1+exp(-HR.95H))
  	  }else{
  	    HR.95L = 0
  	    HR.95H = logb(10)/(1+exp(-HR.95H))
  	  }
  	}else if(HR.95H > 5){
  	  HR = 5/(1+exp(-HR))
  	  if(HR.95L < 0){
  	    if(HR.95L < logb(0.1)){
  	      HR.95L = -0.1/(1+exp(HR.95L))
  	    }
  	    if(HR < 1){
  	      HR.95L = logb(HR) + HR.95L
  	    }  
  	    HR.95H= logb(10)/(1+exp(-HR.95H))
  	  }else{
  	    HR.95L = logb(2)/(1+exp(-HR.95L))
  	    HR.95H = logb(10)/(1+exp(-HR.95H))
  	  }
  	}
  	if(HR.95L < logb(0.1)){
  	  HR.95L = -0.1/(1+exp(HR.95L))
  	  if(HR < 1){
  	    HR.95L =  logb(HR)+HR.95L
  	  }               
  	  if(HR.95H > 10){
  	    HR= 5/(1+exp(-HR))
  	    HR.95H= logb(10)/(1+exp(-HR.95H))
  	  }

  	}
  	if(HR.95H > 3){
  	  HR = 5/(1+exp(-HR))
  	  HR.95H = logb(10)/(1+exp(-HR.95H))
  	  if(HR.95L < 0){
  	    if(HR.95L < logb(0.1)){
  	      HR.95L = -0.1/(1+exp(HR.95L))
  	    }
  	    if(HR < 1){
  	      HR.95L = logb(HR) + HR.95L
  	    }  
  	    HR.95H= logb(10)/(1+exp(-HR.95H))
  	  }else{
  	    HR.95L = logb(2)/(1+exp(-HR.95L))
  	    HR.95H = logb(10)/(1+exp(-HR.95H))
  	  }
  	}
  # 	if(HR<0.001){
  # 	    HR <- 0.005/(1+exp(-HR))
  # 	   
  # 	    if(HR.95H>10){
  # 	    HR.95H <- 10/(1+exp(-HR.95H))
  # 	    }
  # 	    if(HR.95H>5){
  # 	      HR.95H <-  6/(1+exp(-HR.95H))
  # 	    }
  # 	    if(HR.95H>2){
  # 	      HR.95H <-  3/(1+exp(-HR.95H))
  # 	    }
  # 	    if(HR.95H<1){
  # 	      HR.95H <- 0.9/(1+exp(-HR.95H))
  # 	    }
  # 	    if(HR.95L<0){
  # 	    HR.95L <- HR/(1.5+exp(HR.95L))}else{
  # 	      HR.95L <- HR/(1.5+exp(-HR.95L)) 
  # 	      }
  # 	    
  # 	}
  # if(HR.95L>HR){HR.95L=HR/1.5}
  # if(HR.95H<HR){HR.95H=HR*1.5}

  	outTab=rbind(outTab,
  	             cbind(cancer=i,
  	                   HR=HR,
  	                   HR.95L=exp(HR.95L),
  	                   HR.95H=exp(HR.95H),
  	                   pvalue=coxP) )
  }



  # end of loop

  outTab$pvalue=as.numeric(as.character(outTab$pvalue))
  outTab1=arrange(outTab,pvalue)
  significant_cohort=union(significant_cohort,outTab1[which(outTab1$pvalue<0.05),"cancer"])
  write.table(outTab1,file=paste0(gene,".DSS.cox.result.txt"),sep="\t",row.names=F,quote=F)


  bioForest=function(coxFile=null,forestFile=null,forestCol=null){

  	rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  	data=as.data.frame(rt)
  	HR=data[,1:3]
  	hr=sprintf("%.3f",HR[,"HR"])
  	hrLow=sprintf("%.3f",HR[,"HR.95L"])
  	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  	pVal=data[,"pvalue"]
  	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  	clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue") 
  	tabletext <- 
  	  list(c(NA, rownames(HR)),
  	       append("pvalue", pVal),
  	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
  p=forestplot(tabletext, 
               rbind(rep(NA, 3), HR),
               col=clrs,
               graphwidth=unit(50, "mm"),
               xlog=T,
               lwd.ci=1,
               boxsize=0.6,
               xlab="Hazard ratio",
               txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
  )
  	pdf(file=forestFile,width = 9,height = 6,onefile = FALSE)
  	print(p)
  	dev.off()
  }


  bioForest(coxFile=paste0(gene,".DSS.cox.result.txt"),forestFile=paste0(gene,".DSS.forest.pdf"),forestCol="red")
  setwd(ParentWD)
  cp_report(file.path("10.DSS", outdir))
  message("[step_10_cox_dss] done")
}

step_11_cox_dfi <- function() {
  outdir <- paste0("DFI.", gene, ".DFI.univariate_cox_forest_plot")
  message("[step_11_cox_dfi] start")
  setwd(file.path(ParentWD, "11.DFI"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)
  library(limma)
  library(survival)
  library(survminer)
  library(forestplot)

  library(dplyr)
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  files=dir(paste0(ParentWD,"/07.preSurvival/"))
  files=grep(".survival.tsv$",files,value=T)
  files=paste0(ParentWD,"/07.preSurvival/",files)

  surTime=data.frame()
  for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
  }
  colnames(surTime)=c("futime","fustat")
  surTime=surTime[stringr::str_sub(row.names(surTime),-1,-1)=="A",]
  rownames(exp)=exp$samples
  sameSample=intersect(row.names(surTime),row.names(exp))
  surTime=surTime[sameSample,]
  exp=exp[sameSample,]
  surData=cbind(surTime,exp)
  surData=cbind(id=row.names(surData),surData)
  surData=surData[stringr::str_sub(row.names(surData),-3,-2)=="01",]

  library(survival)
  library(survminer)

  surPlot=list()
  rt=surData
  rt$futime=rt$futime/365
  gene=colnames(rt)[7]
  pFilter=0.05
  rt[,"cancer_type"]=factor(rt[,"cancer_type"])



  cliFile="../Survival_SupplementalTable_S1_20171025_xena_sp"

  cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)
  cli=cli[,c("DFI.time","DFI")]
  cli=na.omit(cli)
  colnames(cli)=c("futime","fustat")
  rt$cancer_type=factor(rt$cancer_type)

  outTab=data.frame()
  for(i in levels(rt[,"cancer_type"])){

    rt1=rt[(rt[,"cancer_type"]==i),]
  	data=cbind(rt1,gene=rt1[,gene])
  	data=as.matrix(data[,c(gene,"gene")])
  	if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
  		row.names(data)=gsub(".$","",row.names(data))}
  	data=avereps(data)
  	sameSample=intersect(row.names(data),row.names(cli))
  	if(length(sameSample)<3){next}
  	sameData=data[sameSample,]
  	sameCli=cli[sameSample,]
  	rt1=cbind(sameCli,sameData)
  	rt1$futime=rt1$futime/365


  	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  	pValue=1-pchisq(diff$chisq,df=1)
  	message(paste0("analysisng: ",i," pvalue:",pValue))
  	if(pValue<pFilter){
  		if(pValue<0.001){
  			pValue="p<0.001"
  		}else{
  			pValue=paste0("p=",sprintf("%.03f",pValue))
  		}
  		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)

  		surPlot=ggsurvplot(fit, 
  				    data=rt1,
  				    title=paste0("Cancer: ",i),
  				    pval=pValue,
  				    pval.size=6,
  				    legend.labs=c("high","low"),
  				    legend.title=paste0(gene," levels"),
  				    font.legend=12,
  				    xlab="Time(years)",
  				    ylab="Disease-free interval",
  				    break.time.by = 1,
  				    palette=c("red","blue"),
  				    conf.int=F,
  				    fontsize=4,
  				    risk.table=F)
  		pdf(file=paste0(gene,".DFI.",i,".pdf"),onefile = FALSE, width = 6, height =6)
  		print(surPlot)
  		dev.off()
  	}

  	cox=coxph(Surv(futime, fustat) ~ gene, data = rt1)
  	coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      # export for downstream processing
      HR=coxSummary$conf.int[,"exp(coef)"]
      HR.95L=confint(cox)[,'2.5 %']
      HR.95H=confint(cox)[,'97.5 %']

      # clip values to keep points within the axis range
      if(HR > 5){
        HR = 5/(1+exp(-HR))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else if(HR.95L > 0){
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(7)/(1+exp(-HR.95H))
        }else{
          HR.95L = 0
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }else if(HR.95H > 5){
        HR = 5/(1+exp(-HR))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          if(HR < 1){
            HR.95L = logb(HR) + HR.95L
          }  
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else{
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }
      if(HR.95L < logb(0.1)){
        HR.95L = -0.1/(1+exp(HR.95L))
        if(HR < 1){
          HR.95L =  logb(HR)+HR.95L
        }               
        if(HR.95H > 10){
          HR= 5/(1+exp(-HR))
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }

      }
      if(HR.95H > 3){
        HR = 5/(1+exp(-HR))
        HR.95H = logb(10)/(1+exp(-HR.95H))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          if(HR < 1){
            HR.95L = logb(HR) + HR.95L
          }  
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else{
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }
      outTab=rbind(outTab,
                   cbind(cancer=i,
                         HR=HR,
                         HR.95L=exp(HR.95L),
                         HR.95H=exp(HR.95H),
                         pvalue=coxP) )
  }
  outTab$pvalue=as.numeric(as.character(outTab$pvalue))
  outTab1=arrange(outTab,pvalue)

  significant_cohort=union(significant_cohort,outTab1[which(outTab1$pvalue<0.05),"cancer"])
  write.table(outTab1,file=paste0(gene,".DFI.cox.result.txt"),sep="\t",row.names=F,quote=F)   

  bioForest=function(coxFile=null,forestFile=null,forestCol=null){

  	dff=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  	data=as.data.frame(dff)
  	HR=data[,1:3]
  	hr=sprintf("%.3f",HR[,"HR"])
  	hrLow=sprintf("%.3f",HR[,"HR.95L"])
  	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  	pVal=data[,"pvalue"]
  	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))

  	tabletext <- 
  	  list(c(NA, rownames(HR)),
  	       append("pvalue", pVal),
  	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
  	p=forestplot(tabletext, 
  	             rbind(rep(NA, 3), HR),
  	             # col=clrs,
  	             graphwidth=unit(60, "mm"),
  	             xlog=T,
  	             lwd.ci=1,
  	             boxsize=0.8,
  	             xlab="Hazard ratio",
  	             txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))	)|>fp_set_style(box = "#73038F",
  	                                                                                       line = "black", 
  	                                                                                       summary = "red",align = "lrrr",hrz_lines = "#999999")|>
  	  fp_set_zebra_style("#EFEFEF")
  	pdf(file=forestFile,width = 9,height = 6,onefile = FALSE)
  print(p)
  	dev.off()
  }


  bioForest(coxFile=paste0(gene,".DFI.cox.result.txt"),forestFile=paste0(gene,".DFI.forest.pdf"),forestCol="#ff6f00")
  setwd(ParentWD)
  cp_report(file.path("11.DFI", outdir))
  message("[step_11_cox_dfi] done")
}

step_12_cox_pfi <- function() {
  outdir <- paste0("PFI.", gene, ".PFI.univariate_cox_forest_plot")
  message("[step_12_cox_pfi] start")
  setwd(file.path(ParentWD, "12.PFI"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(limma)
  library(survival)
  library(survminer)
  library(forestplot)

  library(dplyr)
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  files=dir(paste0(ParentWD,"/07.preSurvival/"))
  files=grep(".survival.tsv$",files,value=T)
  files=paste0(ParentWD,"/07.preSurvival/",files)

  surTime=data.frame()
  for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
  }
  colnames(surTime)=c("futime","fustat")
  surTime=surTime[stringr::str_sub(row.names(surTime),-1,-1)=="A",]
  rownames(exp)=exp$samples
  sameSample=intersect(row.names(surTime),row.names(exp))
  surTime=surTime[sameSample,]
  exp=exp[sameSample,]
  surData=cbind(surTime,exp)
  surData=cbind(id=row.names(surData),surData)
  surData=surData[stringr::str_sub(row.names(surData),-3,-2)=="01",]

  library(survival)
  library(survminer)

  surPlot=list()
  rt=surData
  rt$futime=rt$futime/365
  gene=colnames(rt)[7]
  pFilter=0.05
  rt[,"cancer_type"]=factor(rt[,"cancer_type"])

  cliFile="../Survival_SupplementalTable_S1_20171025_xena_sp"

  cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)
  cli=cli[,c("PFI.time","PFI")]
  cli=na.omit(cli)
  colnames(cli)=c("futime","fustat")
  rt$cancer_type=factor(rt$cancer_type)

  outTab=data.frame()
  for(i in levels(rt[,"cancer_type"])){
    rt1=rt[(rt[,"cancer_type"]==i),]

  	data=cbind(rt1,gene=rt1[,gene])
  	data=as.matrix(data[,c(gene,"gene")])
  	if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
  		row.names(data)=gsub(".$","",row.names(data))}
  	data=avereps(data)
  	sameSample=intersect(row.names(data),row.names(cli))
  	sameData=data[sameSample,]
  	sameCli=cli[sameSample,]
  	rt1=cbind(sameCli,sameData)
  	rt1$futime=rt1$futime/365
  	if(nrow(rt1)<3){next}

  	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  	pValue=1-pchisq(diff$chisq,df=1)
  	if(pValue<pFilter){
  		if(pValue<0.001){
  			pValue="p<0.001"
  		}else{
  			pValue=paste0("p=",sprintf("%.03f",pValue))
  		}
  		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)

  		surPlot=ggsurvplot(fit, 
  				    data=rt1,
  				    title=paste0("Cancer: ",i),
  				    pval=pValue,
  				    pval.size=6,
  				    legend.labs=c("high","low"),
  				    legend.title=paste0(gene," levels"),
  				    font.legend=12,
  				    xlab="Time(years)",
  				    ylab="Progression-free interval",
  				    break.time.by = 1,
  				    palette=c("#DBA507","#0D6986"),
  				    conf.int=F,
  				    fontsize=4,
  				    risk.table=FALSE)
  		pdf(file=paste0(gene,".PFI.",i,".pdf"),onefile = FALSE, width = 6, height =6)
  		print(surPlot)
  		dev.off()
  	}

  	cox=coxph(Surv(futime, fustat) ~ gene, data = rt1)
  	coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      # export for downstream processing
      HR=coxSummary$conf.int[,"exp(coef)"]
      HR.95L=confint(cox)[,'2.5 %']
      HR.95H=confint(cox)[,'97.5 %']

      # clip values to keep points within the axis range
      if(HR > 5){
        HR = 5/(1+exp(-HR))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else if(HR.95L > 0){
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(7)/(1+exp(-HR.95H))
        }else{
          HR.95L = 0
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }else if(HR.95H > 5){
        HR = 5/(1+exp(-HR))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          if(HR < 1){
            HR.95L = logb(HR) + HR.95L
          }  
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else{
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }
      if(HR.95L < logb(0.1)){
        HR.95L = -0.1/(1+exp(HR.95L))
        if(HR < 1){
          HR.95L =  logb(HR)+HR.95L
        }               
        if(HR.95H > 10){
          HR= 5/(1+exp(-HR))
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }

      }
      if(HR.95H > 3){
        HR = 5/(1+exp(-HR))
        HR.95H = logb(10)/(1+exp(-HR.95H))
        if(HR.95L < 0){
          if(HR.95L < logb(0.1)){
            HR.95L = -0.1/(1+exp(HR.95L))
          }
          if(HR < 1){
            HR.95L = logb(HR) + HR.95L
          }  
          HR.95H= logb(10)/(1+exp(-HR.95H))
        }else{
          HR.95L = logb(2)/(1+exp(-HR.95L))
          HR.95H = logb(10)/(1+exp(-HR.95H))
        }
      }
      outTab=rbind(outTab,
                   cbind(cancer=i,
                         HR=HR,
                         HR.95L=exp(HR.95L),
                         HR.95H=exp(HR.95H),
                         pvalue=coxP) )
  }

  outTab$pvalue=as.numeric(as.character(outTab$pvalue))
  outTab1=arrange(outTab,pvalue)
  significant_cohort=union(significant_cohort,outTab1[which(outTab1$pvalue<0.05),"cancer"])
  write.table(outTab1,file=paste0(gene,".PFI.cox.result.txt"),sep="\t",row.names=F,quote=F)   


  bioForest=function(coxFile=null,forestFile=null,forestCol=null){

  	rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  	data=as.data.frame(rt)
  	HR=data[,1:3]
  	hr=sprintf("%.3f",HR[,"HR"])
  	hrLow=sprintf("%.3f",HR[,"HR.95L"])
  	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  	pVal=data[,"pvalue"]
  	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  	clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue")     
  	tabletext <- 
  	  list(c(NA, rownames(HR)),
  	       append("pvalue", pVal),
  	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) ) 
   p=forestplot(tabletext, 
                rbind(rep(NA, 3), HR),
                col=clrs,
                graphwidth=unit(50, "mm"),
                xlog=T,
                lwd.ci=1,
                boxsize=0.6,
                xlab="Hazard ratio",
                txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25)))
  	pdf(file=forestFile,width = 9,height = 6,onefile = FALSE)
  print(p)
  	dev.off()
  }


  bioForest(coxFile=paste0(gene,".PFI.cox.result.txt"),forestFile=paste0(gene,".PFI.forest.pdf"),forestCol="#2f9e44")
  setwd(ParentWD)
  cp_report(file.path("12.PFI", outdir))
  message("[step_12_cox_pfi] done")
}

step_16_tmb_radar <- function() {
  outdir <- paste0("TMBradar.", gene, ".TMB.radarplot")
  message("[step_16_tmb_radar] start")
  setwd(file.path(ParentWD, "16.TMBradar"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)


  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  rownames(exp)=exp$samples
  TMB=read.table(paste0(ParentWD,"/15.TMBcor/TMB.txt"), header=T,sep="\t",row.names=1,check.names=F)
  exp$cancer_type=factor(exp$cancer_type)

  group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  exp=exp[group==0,]

  sameSample=intersect(row.names(TMB),row.names(exp))
  TMB=TMB[sameSample,]
  exp=exp[sameSample,]


  outTab=data.frame()
  fmsbTab=data.frame()
  for(i in levels(exp[,"cancer_type"])){
    exp1=exp[(exp[,"cancer_type"]==i),]
    TMB1=TMB[(TMB[,"CancerType"]==i),]
    x=as.numeric(TMB1[,1])
    y=as.numeric(exp1[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
    outTab=rbind(outTab,cbind(CancerType=i,cor=cor,pValue=pValue,sig))
    fmsbTab=rbind(fmsbTab,cbind(CancerType=i,cor=cor))
  }
  write.table(outTab,file=paste0(gene,".TMB.radar.corStat.txt"),sep="\t",row.names=F,quote=F)
  write.table(t(fmsbTab),file=paste0(gene,".TMB.radar.fmsbInput.txt"),sep="\t",col.names=F,quote=F)


  library(fmsb) 

  data=read.table(paste0(gene,".TMB.radar.fmsbInput.txt"),sep = "\t",header = T,row.names = 1)
  maxValue=ceiling(max(abs(data))*10)/10
  data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)

  colors="#ff5722"

  corStat=read.table(paste0(gene,".TMB.radar.corStat.txt"),header=T,sep="\t",row.names=1,check.names=F)
  colnames(data)=paste0(colnames(data),corStat$sig)


  pdf(file=paste0(gene,".TMB.cor.radar.pdf"),height=7,width=7)
  radarchart( data, axistype=1 , 
      pcol=colors,                
      plwd=2 ,                  
      plty=1,                   
      cglcol="#bcaaa4",             
      cglty=1,                  
      caxislabels=seq(-maxValue,maxValue,maxValue/2),   
      cglwd=1.2,                 
      axislabcol="#1e88e5",  
      vlcex=0.8            
  )
  dev.off()
  setwd(ParentWD)
  cp_report(file.path("16.TMBradar", outdir))
  message("[step_16_tmb_radar] done")
}

step_18_msi_radar <- function() {
  outdir <- paste0("MSIradar.", gene, ".MSI.radarplot")
  message("[step_18_msi_radar] start")
  setwd(file.path(ParentWD, "18.MSIradar"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  rownames(exp)=exp$samples
  exp$cancer_type=factor(exp$cancer_type)
  MSI=read.table(paste0(ParentWD,"/17.MSIcor/MSI.txt"), header=T,sep="\t",row.names=1,check.names=F)

  group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  exp=exp[group==0,]

  sameSample=intersect(row.names(MSI),row.names(exp))
  MSI=MSI[sameSample,]
  exp=exp[sameSample,]

  exp$cancer_type=factor(exp$cancer_type)
  outTab=data.frame()
  fmsbTab=data.frame()
  for(i in levels(exp[,"cancer_type"])){
    exp1=exp[(exp[,"cancer_type"]==i),]
    MSI1=MSI[(MSI[,"CancerType"]==i),]
    x=as.numeric(MSI1[,1])
    y=as.numeric(exp1[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
    outTab=rbind(outTab,cbind(CancerType=i,cor=cor,pValue=pValue,sig))
    fmsbTab=rbind(fmsbTab,cbind(CancerType=i,cor=cor))
  }
  write.table(outTab,file=paste0(gene,".msi.corStat.txt"),sep="\t",row.names=F,quote=F)
  write.table(t(fmsbTab),file=paste0(gene,".msi.fmsbInput.txt"),sep="\t",col.names=F,quote=F)


  library(fmsb) 

  data=read.table(paste0(gene,".msi.fmsbInput.txt"),header=T,sep="\t",row.names=1,check.names=F)  
  maxValue=ceiling(max(abs(data))*10)/10
  data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)

  colors="blue"

  corStat=read.table(paste0(gene,".msi.corStat.txt"),header=T,sep="\t",row.names=1,check.names=F)
  colnames(data)=paste0(colnames(data),corStat$sig)


  pdf(file=paste0(gene,".MSI.cor.radar.pdf"),height=7,width=7)
  radarchart( data, axistype=1 , 
      pcol=colors,                
      plwd=2 ,                     
      plty=1,                      
      cglcol="grey",               
      cglty=1,                      
      caxislabels=seq(-maxValue,maxValue,maxValue/2), 
      cglwd=1.2,                  
      axislabcol="green",       
      vlcex=0.8       
  )
  dev.off()
  setwd(ParentWD)
  cp_report(file.path("18.MSIradar", outdir))
  message("[step_18_msi_radar] done")
}

step_20_estimate <- function() {
  outdir <- paste0("9.", gene, ".Estimate.Score.corplot")
  message("[step_20_estimate] start")
  setwd(file.path(ParentWD, "20.estimateCor"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(ggplot2)
  library(ggpubr)
  library(ggExtra)

  pFilter=0.01  

  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"

  rownames(exp)=exp$samples
  exp$cancer_type=factor(exp$cancer_type)
  TME=read.table(paste0(ParentWD,"/20.estimateCor/estimateScores.txt"), header=T,sep="\t",row.names=1,check.names=F)

  group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  exp=exp[group==0,]


  TME=TME[stringr::str_sub(row.names(TME),-1,-1)=="A",]
  sameSample=intersect(row.names(TME),row.names(exp))
  TME=TME[sameSample,]
  exp=exp[sameSample,]


  outTab=data.frame()

  for(i in levels(exp[,"cancer_type"])){
      exp1=exp[(exp[,"cancer_type"]==i),]
      TME1=TME[(TME[,"CancerType"]==i),]
      y=as.numeric(exp1[,gene])
      outVector=data.frame(i,gene)

  	for(j in colnames(TME1)[1:2]){
  		x=as.numeric(TME1[,j])
  		df1=as.data.frame(cbind(x,y))
  		corT=cor.test(x,y,method="spearman")
  		cor=corT$estimate
  		pValue=corT$p.value
  		outVector=cbind(outVector,pValue)
  		p1=ggplot(df1, aes(x, y)) + 
  			xlab(j)+ylab(gene)+
  			ggtitle(paste0("Cancer: ",i))+theme(title=element_text(size=10))+
  		    geom_point()+ geom_smooth(method="lm") + theme_bw()+
  		    stat_cor(method = 'spearman', aes(x =x, y =y))
  	    p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
  		if(pValue<pFilter){
  			pdf(file=paste0(gene,".estimateCor.",i,"_",j,".pdf"),width=5,height=5)
  			print(p2)
  			dev.off()
  		}
  	}
  	outTab=rbind(outTab,outVector)
  }
  colNames=c("CancerType","Gene",colnames(TME)[1:2])
  colnames(outTab)=colNames
  write.table(outTab,file=paste0(gene,".estimateCor.result.txt"),sep="\t",row.names=F,quote=F)
  setwd(ParentWD)
  cp_report(file.path("20.estimateCor", outdir))
  message("[step_20_estimate] done")
}

step_21_cibersort <- function() {
  outdir <- paste0("x10.", gene, ".Cibersort.statplots")
  message("[step_21_cibersort] start")
  setwd(file.path(ParentWD, "21.CIBERSORT"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)


  library(ggpubr)

  pFilter=0.05
  TME=read.table(paste0(ParentWD,"/21.CIBERSORT/CIBERSORT.result.txt"),sep="\t",header=T,row.names=1,check.names=F)
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"
  rownames(exp)=exp$samples
  exp$cancer_type=factor(exp$cancer_type)
  group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  exp=exp[group==0,]
  sameSample=intersect(row.names(TME),row.names(exp))
  TME=TME[sameSample,]
  exp=exp[sameSample,]
  TME$CancerType=factor(TME$CancerType)
  TME[,"gene"]=exp[,gene]
  outTab=data.frame(logFCLH=rep(1,ncol(TME)-2),logFCHL=rep(1,ncol(TME)-2),pvalue=rep(1,ncol(TME)-2),row.names = colnames(TME)[1:(ncol(TME)-2)])
  for(i in levels(TME$CancerType)){

  data=TME[(TME[,"CancerType"]==i),]
  Low=data[which(data$gene<=median(data$gene)),]
  High=data[which(data$gene>median(data$gene)),]
  High=High%>%dplyr::select(-gene,-CancerType)
  Low=Low%>%dplyr::select(-gene,-CancerType)
  n=1
  for(j in rownames(outTab)){
    outTab$logFCLH[n]=log((mean(as.numeric(Low[,j]),na.rm=T)+0.0001)/(mean(as.numeric(High[,j]),na.rm=T)+0.0001),2)
    outTab$logFCHL[n]=log((mean(as.numeric(High[,j]),na.rm=T)+0.0001)/(mean(as.numeric(Low[,j]),na.rm=T)+0.0001),2)
    outTab$pvalue[n]=wilcox.test(as.numeric(Low[,j]),as.numeric(High[,j]))$p.value
    n=n+1
  }


  outTab$sig=ifelse(outTab$pvalue<0.05,"yes","no")
  outTab$celltype=rownames(outTab)
  p=ggbarplot(outTab,x="celltype",y="logFCLH",fill = "sig", palette = c("#d9d9d9", "#b71c1c"),xlab = "",ylab = "LogFC")+rotate_x_text(60)

  pdf(file=paste0(gene,"_LowHigh_",i,"_barplot_cibersort.pdf"),width=6,height=4)
  print(p)
  dev.off()

  p1=ggbarplot(outTab,x="celltype",y="logFCHL",fill = "sig", palette = c("#d9d9d9", "#b71c1c"),xlab = "",ylab = "LogFC")+rotate_x_text(60)
  pdf(file=paste0(gene,"_HighLow_",i,"_barplot_cibersort.pdf"),width=6,height=4)
  print(p1)
  dev.off()

  High=t(as.matrix(High))
  Low=t(as.matrix(Low))
  col1=rainbow(nrow(High),s=0.7,v=0.7)
  col2=rainbow(nrow(Low),s=0.7,v=0.7)
  pdf(paste0(gene,".",i,".High.cibersort.pdf"),height=10,width=22)
  par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
  a1 = barplot(High,col=col1,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
  a2=axis(2,tick=F,labels=F)
  axis(2,a2,paste0(a2*100,"%"))
  axis(1,a1,labels=F)
  par(srt=60,xpd=T);text(a1,-0.02,colnames(High),adj=1,cex=0.6);par(srt=0)
  ytick2 = cumsum(High[,ncol(High)])
  ytick1 = c(0,ytick2[-length(ytick2)])
  legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(High),col=col1,pch=15,bty="n",cex=1.3)
  dev.off()
  pdf(paste0(gene,".",i,".Low.cibersort.pdf"),height=10,width=22)
  par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
  a1 = barplot(Low,col=col1,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
  a2=axis(2,tick=F,labels=F)
  axis(2,a2,paste0(a2*100,"%"))
  axis(1,a1,labels=F)
  par(srt=60,xpd=T);text(a1,-0.02,colnames(Low),adj=1,cex=0.6);par(srt=0)
  ytick2 = cumsum(Low[,ncol(Low)])
  ytick1 = c(0,ytick2[-length(ytick2)])
  legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(Low),col=col2,pch=15,bty="n",cex=1.3)
  dev.off()

  }
  setwd(ParentWD)
  cp_report(file.path("21.CIBERSORT", outdir))
  message("[step_21_cibersort] done")
}

step_22_cibersort_cor <- function() {
  outdir <- paste0("x11.", gene, ".Cibersort.Corr")
  message("[step_22_cibersort_cor] start")
  setwd(file.path(ParentWD, "22.CIBERSORTcor"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(ggplot2)
  library(ggpubr)
  library(ggExtra)

  pFilter=0.01  
  TME=read.table(paste0(ParentWD,"/21.CIBERSORT/CIBERSORT.result.txt"),sep="\t",header=T,row.names=1,check.names=F)
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T)
  colnames(exp)[1]="samples"
  rownames(exp)=exp$samples
  exp$cancer_type=factor(exp$cancer_type)
  group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  exp=exp[group==0,]
  sameSample=intersect(row.names(TME),row.names(exp))
  TME=TME[sameSample,]
  exp=exp[sameSample,]


  outTab=data.frame()

  for(i in levels(exp[,"cancer_type"])){
      exp1=exp[(exp[,"cancer_type"]==i),]
      TME1=TME[(TME[,"CancerType"]==i),]
      y=as.numeric(exp1[,gene])
      outVector=data.frame(i,gene)

  	for(j in colnames(TME1)[1:22]){
  		x=as.numeric(TME1[,j])
  		if(sd(x)>0.01){
  			df1=as.data.frame(cbind(x,y))
  			corT=cor.test(x,y,method="spearman")
  			cor=corT$estimate
  			pValue=corT$p.value
  			outVector=cbind(outVector,pValue)
  			p1=ggplot(df1, aes(x, y)) + 
  				xlab(j)+ylab(gene)+
  				ggtitle(paste0("Cancer: ",i))+theme(title=element_text(size=10))+
  			    geom_point()+ geom_smooth(method="lm") + theme_bw()+
  			    stat_cor(method = 'spearman', aes(x =x, y =y))
  		    p2=ggMarginal(p1, type = "density", xparams = list(fill = "#a4b0be"),yparams = list(fill = "#d35400"))
  			if(pValue<pFilter){
  				pdf(file=paste0(gene,".CIBERSORTCor.",i,"_",j,".pdf"),width=5,height=5)
  				print(p2)
  				dev.off()
  			}
  		}
  		else{
  			outVector=cbind(outVector,pValue=1)
  		}
  	}
  	outTab=rbind(outTab,outVector)
  }
  colNames=c("CancerType","Gene",colnames(TME)[1:22])
  colnames(outTab)=colNames
  write.table(outTab,file=paste0(gene,".CIBERSORTcor.result.txt"),sep="\t",row.names=F,quote=F)
  setwd(ParentWD)
  cp_report(file.path("22.CIBERSORTcor", outdir))
  message("[step_22_cibersort_cor] done")
}

step_26_immunefunction <- function() {
  outdir <- paste0("x26.", gene, ".immunefunction")
  message("[step_26_immunefunction] start")
  setwd(file.path(ParentWD, "26.immunefunction"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  #grouplist: sample*group
  #x: Cohort  sample*gene(or risk)
  #y: Cohort  sample*multiple GSVA
  library(limma)
  library(ggpubr)
  ScoreFile="immune.GSVA.pancancer.txt"
  ScoreType=gsub("\\.GSVA.*txt$","",ScoreFile)# analysis type label
  Score=read.table(paste0(ParentWD,"/database/",ScoreFile),sep="\t",header = T,row.names = 1)      # input
  #Patient  Cohort
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T,row.names = 1)
  PaintScore=gene

  GSVARelationship=function(x,y,Tumors="ALL"){
  if(Tumors=="ALL"){
  Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  }
  valueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
  rownames(valueDF)=colnames(y)[3:ncol(y)]
  colnames(valueDF)=Ts
  pvalueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
  rownames(pvalueDF)=colnames(y)[3:(ncol(y))]
  colnames(pvalueDF)=Ts
  for(j in Ts){
  if(!(j %in% y$Cohort) | !(j %in% x$cancer_type)){
  message(paste0("Seems ",j," is not present in neither of input data"))
  next
  }

  x1=x[which(x$cancer_type==j),]
  y1=y[which(y$Cohort==j),]
  intersample=intersect(rownames(x1),rownames(y1))
  x1=x1[intersample,]
  y1=y1[intersample,]
  if(length(intersample)<nrow(x1)/3){
  message(paste0("Seems ",i," metrics of ", j," is not enough")) 
  }
  for(i in colnames(y)[2:ncol(y)]){
  tmp=cor.test(x1[,PaintScore],y1[,i])
  valueDF[i,j]=tmp$estimate
  pvalueDF[i,j]=tmp$p.value
  }
  message(paste0("Done ",j))
  }
  list_data <- list(valueDF,pvalueDF)
  names(list_data) <- c("value","pvalue")
  return(list_data)
  } 

  heatmapData=GSVARelationship(exp,Score)# x: PaintScore (expression), y: pathway/cell-type scores
  hmpvalue <- as.data.frame(heatmapData$pvalue) 
  hmdata <- as.data.frame(heatmapData$value)  

  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        hmdata[i,j]=NA
      }
    }
  }
  hmdata=zz(hmdata)
  hmpvalue=zz(hmpvalue)
  Cs.mat=data.frame(matrix(nrow = nrow(hmpvalue),ncol=ncol(hmpvalue)))
  rownames(Cs.mat)=row.names(hmpvalue)
  colnames(Cs.mat)=colnames(hmpvalue)
  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        Cs.mat[i,j]='X'
      }else{Cs.mat[i,j]=''}
    }
  }

  corl=c("#244B9E", "white","firebrick")
  a=pheatmap(hmdata,cluster_rows = F,cluster_cols = F,na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
  pdf(paste0(gene,".",ScoreType,".scoreCor.pheatmap.pdf"),width = 12,height = 15)
  print(a)
  dev.off()


  library(reshape2)
  GSVADeplot=function(x,y,Tumors="ALL"){
  if(Tumors=="ALL"){
    Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  }
    plotlist = list()
    count=1
  for(j in Ts){
      if(!(j %in% y$Cohort) | !(j %in% x$cancer_type)){
        message(paste0("Seems ",j," is not present in neither of input data"))
        next
      }

      x1=x[which(x$cancer_type==j),]
      y1=y[which(y$Cohort==j),]
      intersample=intersect(rownames(x1),rownames(y1))
      y1=y1[intersample,]
      if(length(intersample)<nrow(x1)/3){
        message(paste0("Seems ",i," metrics of ", j," is not enough")) 
      }

  x1=x1[intersample, PaintScore, drop=F]
  data=cbind(x1, y1)
  data=data%>%dplyr::select(-Cohort)
  colnames(data)[1]="Risk"
  data$Risk=ifelse(data$Risk>median(data$Risk),"High", "Low")
  group=levels(factor(data$Risk))
  data$Risk=factor(data$Risk, levels=c("Low", "High"))
  group=levels(factor(data$Risk))
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  for(i in 2:ncol(data)){data[,i]=as.numeric(data[,i])}
  rt=reshape2::melt(data,id.vars="Risk")
  colnames(rt)=c(PaintScore,"Type","Score")
  p <- ggboxplot(rt, x="Type", y="Score", color = PaintScore,main=j,ylab="Score",add = "none",xlab="",palette = c("#2762ae","#e71b1c") )+rotate_x_text(50)+stat_compare_means(aes(group=rt[,PaintScore]),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  ggsave(paste0(PaintScore,".group.",ScoreType,".",j,".pdf"), plot = p,width=5+0.3*ncol(data), height=5)
  count=count+1
  }
  #pdf(file=paste0(PaintScore,".group.",ScoreType,".",j,".pdf"), width=5+0.3*ncol(data), height=5)
  #p
  #invisible(dev.off())

  }
  GSVADeplot(exp,Score)
  setwd(ParentWD)
  cp_report(file.path("26.immunefunction", outdir))
  message("[step_26_immunefunction] done")
}

step_27_timer2 <- function() {
  outdir <- paste0("x27.", gene, ".timer2infiltration")
  message("[step_27_timer2] start")
  setwd(file.path(ParentWD, "27.infiltration"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  #grouplist: sample*group
  #x: Cohort  sample*gene(or risk)
  #y: Cohort  sample*multiple GSVA
  library(limma)
  library(ggpubr)
  GSVARelationship=function(x,y,Tumors="ALL"){
    if(Tumors=="ALL"){
      Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
    }
    valueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
    rownames(valueDF)=colnames(y)[3:ncol(y)]
    colnames(valueDF)=Ts
    pvalueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
    rownames(pvalueDF)=colnames(y)[3:(ncol(y))]
    colnames(pvalueDF)=Ts
    for(j in Ts){
      if(!(j %in% y$Cohort) | !(j %in% x$cancer_type)){
        message(paste0("Seems ",j," is not present in neither of input data"))
        next
      }

      x1=x[which(x$cancer_type==j),]
      y1=y[which(y$Cohort==j),]
      intersample=intersect(rownames(x1),rownames(y1))
      x1=x1[intersample,]
      y1=y1[intersample,]
      if(length(intersample)<nrow(x1)/3){
        message(paste0("Seems ",i," metrics of ", j," is not enough")) 
      }
      for(i in colnames(y)[2:ncol(y)]){
        if(sum(is.na(y1[,i]))==nrow(y1)){next}
        tmp=cor.test(x1[,PaintScore],y1[,i])
        valueDF[i,j]=tmp$estimate
        pvalueDF[i,j]=tmp$p.value
      }
      message(paste0("Done ",j))
    }
    list_data <- list(valueDF,pvalueDF)
    names(list_data) <- c("value","pvalue")
    return(list_data)
  } 
  library(gg.gap)
  library(reshape2)
  GSVADeplot=function(x,y,Tumors="ALL"){
    if(Tumors=="ALL"){
      Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
    }
    plotlist = list()
    count=1
    for(j in Ts){
      if(!(j %in% y$Cohort) | !(j %in% x$cancer_type)){
        message(paste0("Seems ",j," is not present in neither of input data"))
        next
      }

      x1=x[which(x$cancer_type==j),]
      y1=y[which(y$Cohort==j),]
      intersample=intersect(rownames(x1),rownames(y1))
      y1=y1[intersample,]
      if(length(intersample)<nrow(x1)/3){
        message(paste0("Seems ",i," metrics of ", j," is not enough")) 
      }

      x1=x1[intersample, PaintScore, drop=F]
      data=cbind(x1, y1)
      data=data%>%dplyr::select(-Cohort)
      colnames(data)[1]="Risk"
      data$Risk=ifelse(data$Risk>median(data$Risk),"High", "Low")
      group=levels(factor(data$Risk))
      data$Risk=factor(data$Risk, levels=c("Low", "High"))
      group=levels(factor(data$Risk))
      comp=combn(group,2)
      my_comparisons=list()
      for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
      for(i in 2:ncol(data)){data[,i]=as.numeric(data[,i])}
      rt=reshape2::melt(data,id.vars="Risk")
      colnames(rt)=c(PaintScore,"Type","Score")
      trans <- function(x){pmin(x,3) + 0.005*pmax(x-3,0)}
      rt$Score1=trans(rt$Score)
      p <- ggboxplot(rt, x="Type", y="Score1", color = PaintScore,main=j,ylab="Score",add = "none",xlab="",palette = c("#2762ae","#e71b1c") )+rotate_x_text(50)+coord_cartesian(ylim = c(0,7))+stat_compare_means(aes(group=rt[,PaintScore]),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",label.y = 7)
      ggsave(paste0(PaintScore,".group.",ScoreType,".",j,".pdf"), plot = p,width=5+0.3*ncol(data), height=6)
      count=count+1
    }
    #pdf(file=paste0(PaintScore,".group.",ScoreType,".",j,".pdf"), width=5+0.3*ncol(data), height=5)
    #p
    #invisible(dev.off())

  }
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T,row.names = 1)
  PaintScore=gene
  myfiles=c("pancancer_Bcell_timer2.txt","pancancer_Tcellcd4_timer2.txt","pancancer_Tcellcd8_timer2.txt","pancancer_Tcellreg_timer2.txt","pancancer_Macrophage_timer2.txt","pancancer_fibroblast_timer2.txt","pancancer_score_timer2.txt")
  for( mf in myfiles){
  ScoreFile=mf
  ScoreType=stringr::str_match(ScoreFile,"pancancer_(.*)_timer2.*")[2]# analysis type label
  dir.create(ScoreType)
  setwd(ScoreType)
  Score=read.table(paste0(ParentWD,"/database/",ScoreFile),sep="\t",header = T,row.names = 1)      # input
  #Patient  Cohort
  heatmapData=GSVARelationship(exp,Score)# x: PaintScore (expression), y: pathway/cell-type scores
  hmpvalue <- as.data.frame(heatmapData$pvalue) 
  hmdata <- as.data.frame(heatmapData$value)  

  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(is.na(hmpvalue[i,j])){
        hmdata[i,j]=NA
        next
      }
      if(hmpvalue[i,j]>=0.05){
        hmdata[i,j]=NA
      }
    }
  }
  hmdata=zz(hmdata)
  hmpvalue=zz(hmpvalue)
  Cs.mat=data.frame(matrix(nrow = nrow(hmpvalue),ncol=ncol(hmpvalue)))
  rownames(Cs.mat)=row.names(hmpvalue)
  colnames(Cs.mat)=colnames(hmpvalue)
  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(is.na(hmpvalue[i,j])){
        Cs.mat[i,j]='X'
        next
      }
      if(hmpvalue[i,j]>=0.05){
        Cs.mat[i,j]='X'
      }else{Cs.mat[i,j]=''}
    }
  }

  corl=c("#244B9E", "white","firebrick")
  a=pheatmap(hmdata,cluster_rows = F,cluster_cols = F,na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
  pdf(paste0(gene,".",ScoreType,".scoreCor.pheatmap.pdf"),width = 12,height = 15)
  print(a)
  dev.off()
  GSVADeplot(exp,Score)
  setwd("..")
  }
  setwd(ParentWD)
  cp_report(file.path("27.infiltration", outdir))
  message("[step_27_timer2] done")
}

step_28_hallmark_gsva <- function() {
  outdir <- paste0("x28.", gene, ".HallmarkGSVA")
  message("[step_28_hallmark_gsva] start")
  setwd(file.path(ParentWD, "28.HallmarkGSVA"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  #grouplist: sample*group
  #x: Cohort  sample*gene(or risk)
  #y: Cohort  sample*multiple GSVA
  library(limma)
  library(ggpubr)


  library(reshape2)
  library(RColorBrewer)
  options(stringsAsFactors = F)

  ScoreFile="hallmark.GSVA.pancancer.txt"
  ScoreType=gsub("\\.GSVA.*txt$","",ScoreFile)# analysis type label
  Score=read.table(paste0(ParentWD,"/database/",ScoreFile),sep="\t",header = T,row.names = 1)      # input
  #Patient  Cohort
  exp=read.table(paste0(ParentWD,"/data/",gene,"exprTime.txt"),sep=",",header = T,row.names = 1)
  PaintScore=gene
  Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  GSVARelationship=function(x,y,Tumors="ALL"){
  if(Tumors=="ALL"){
    Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")#
  }
  valueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
  rownames(valueDF)=colnames(y)[3:ncol(y)]
  colnames(valueDF)=Ts
  pvalueDF=data.frame(matrix(nrow=ncol(y)-2,ncol=length(Ts)))
  rownames(pvalueDF)=colnames(y)[3:(ncol(y))]
  colnames(pvalueDF)=Ts
  for(j in Ts){
  if(!(j %in% y$Cohort) | !(j %in% x$cancer_type)){
  message(paste0("Seems ",j," is not present in neither of input data"))
  next
  }

  x1=x[which(x$cancer_type==j),]
  y1=y[which(y$Cohort==j),]
  intersample=intersect(rownames(x1),rownames(y1))
  x1=x1[intersample,]
  y1=y1[intersample,]
  if(length(intersample)<nrow(x1)/3){
  message(paste0("Seems ",i," metrics of ", j," is not enough")) 
  }
  for(i in colnames(y)[2:ncol(y)]){
  tmp=cor.test(x1[,PaintScore],y1[,i])
  valueDF[i,j]=tmp$estimate
  pvalueDF[i,j]=tmp$p.value
  }
  message(paste0("Done ",j))
  }
  list_data <- list(valueDF,pvalueDF)
  names(list_data) <- c("value","pvalue")
  return(list_data)
  } 

  heatmapData=GSVARelationship(exp,Score,Tumors="ALL")# x: PaintScore (expression), y: pathway/cell-type scores
  hmpvalue <- as.data.frame(heatmapData$pvalue) 
  hmdata <- as.data.frame(heatmapData$value)  

  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        hmdata[i,j]=0
       # hmdata[i,j]=NA
      }
    }
  }

  Cs.mat=data.frame(matrix(nrow = nrow(hmpvalue),ncol=ncol(hmpvalue)))
  rownames(Cs.mat)=row.names(hmpvalue)
  colnames(Cs.mat)=colnames(hmpvalue)
  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        Cs.mat[i,j]='X'
      }else{Cs.mat[i,j]=''}
    }
  }


  Ts_factor=factor(Ts)
  significant_cohort_factor=factor(significant_cohort)
  nosig_cohort_factor=factor(setdiff(Ts,significant_cohort))
  hmdata1=hmdata[,c(levels(significant_cohort_factor),levels(nosig_cohort_factor))]
  corl=c("#51A3CC", "white","#A55D9C")
  a=pheatmap(hmdata1,cluster_rows = T,cluster_cols = F,gaps_col =length(significant_cohort_factor),na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
  pdf(paste0(gene,".",ScoreType,".scoreCor.pheatmap.pdf"),width = 12,height = 12)
  print(a)
  dev.off()
  #a=sort(rowSums(!is.na(hmdata)),decreasing = T)
  #a1=a[a>median(a)]

  a=sort(rowSums(hmdata==0),decreasing = T)
  a1=a[a>median(a)]
  # a1: genes significant in more pathways than the median

  gmtFile="h.all.v7.4.symbols.gmt"  # GMT file
  geneSet=GSEABase::getGmt(paste0(ParentWD,"/database/",gmtFile))
  geneSet1=geneSet[paste0("HALLMARK_",names(a1))]
  for(gs in 1:length(geneSet1)){


  geneRTtype=geneSet1[gs]@.Data[[1]]@setName
  geneRT= geneSet1[gs]@.Data[[1]]@geneIds
  message(paste0("Analysising the following geneset: ",geneRTtype))
  files=dir(paste0(ParentWD,"/data"))
  files=grep("^symbol.*\\.tsv$",files,value=T,perl = T)


  outTab=data.frame()
  corTab=data.frame()
  sameGenes=c()

  pb=txtProgressBar(style=3)# progress bar
  start_time=Sys.time()
  for(i in 1:length(files)){

    CancerType=gsub("^symbol\\.|\\.tsv$","",files[i])
    message(paste0("Analysising cancer type: ",CancerType))
    setTxtProgressBar(pb,i/length(files))
    rt=readRDS(paste0(ParentWD,"/data/",files[i],".RDS"))
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    data=avereps(data)


    group=sapply(strsplit(colnames(data),"\\-"),"[",4)
    group=sapply(strsplit(group,""),"[",1)
    group=gsub("2","1",group)
    data=data[,group==0]

    sameGenes=intersect(geneRT,row.names(data))
    data=data[c(sameGenes,gene),]


    x=as.numeric(data[gene,])
    outVector=data.frame(CancerType)
    corVector=data.frame(CancerType)
    for(j in sameGenes){
      y=as.numeric(data[j,])
      if(sd(y)>0.01){
        corT=cor.test(x,y)
        cor=corT$estimate
        pValue=corT$p.value
        outVector=cbind(outVector,pValue)
        corVector=cbind(corVector,cor)
      }
      else{
        outVector=cbind(outVector,pValue=1)
        corVector=cbind(corVector,cor=0)
      }
    }
    outTab=rbind(outTab,outVector)
    corTab=rbind(corTab,corVector)
  }
  end_time=Sys.time()
  close(pb)
  run_time=end_time-start_time
  # close progress bar

  colnames(outTab)=c("CancerType",sameGenes)
  outTab$CancerType=as.character(outTab$CancerType)
  outTab$CancerType=stringr::str_remove(outTab$CancerType,"\\.htseq_fpkm")

  write.table(outTab,file=paste0(gene,".",geneRTtype,".geneCor.pvalue.txt"),sep="\t",row.names=F,quote=F)


  colnames(corTab)=c("CancerType",sameGenes)
  corTab$CancerType=as.character(corTab$CancerType)
  corTab$CancerType=stringr::str_remove(corTab$CancerType,"\\.htseq_fpkm")
  write.table(corTab,file=paste0(gene,".",geneRTtype,".geneCor.cor.txt"),sep="\t",row.names=F,quote=F)


  hmpvalue <- rn(outTab)
  phmdata <- rn(corTab)


  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        phmdata[i,j]=0
      }
    }
  }
  phmdata=zz(phmdata)
  hmpvalue=zz(hmpvalue)

  if("FAS"%in%rownames(phmdata)){
    phmdata=phmdata[which(rownames(phmdata)!="FAS"),]
    hmpvalue=hmpvalue[which(rownames(hmpvalue)!="FAS"),]
  }



  #a=rownames(phmdata)
  #phmdata=phmdata[sort(a),]

  phmdata_scaled<-t(scale(t(as.matrix(phmdata))))
  phmdata_scaled_dist=dist(phmdata_scaled)
  phmdata_scaled_dist[is.na(phmdata_scaled_dist)] <- 1
  phm_hcluster<-hclust(phmdata_scaled_dist,method = "average")# gene order from average-linkage clustering, used for pagination
  plot(phm_hcluster,hang=-1,cex=.8,main="title")

  phmdata1=phmdata[,c(levels(significant_cohort_factor),levels(nosig_cohort_factor))]
  phmdata1=phmdata1[phm_hcluster$labels[phm_hcluster$order],]
  hmpvalue1=hmpvalue[,c(levels(significant_cohort_factor),levels(nosig_cohort_factor))]
  hmpvalue1=hmpvalue1[phm_hcluster$labels[phm_hcluster$order],]
  Cs.mat=data.frame(matrix(nrow = nrow(hmpvalue1),ncol=ncol(hmpvalue1)))
  rownames(Cs.mat)=row.names(hmpvalue1)
  colnames(Cs.mat)=colnames(hmpvalue1)
  for(i in 1:nrow(hmpvalue1)){
    for(j in 1:ncol(hmpvalue1)){
      if(hmpvalue1[i,j]>=0.05){
        Cs.mat[i,j]='X'
      }else{Cs.mat[i,j]=''}
    }
  }

  corl=c("#396b56", "white","#9d08aa")
  nrow_to_print=50
  n_pages=seq(1,nrow(phmdata1),nrow_to_print)
  for(print_page in 1:length(n_pages)){
    tmp1=(print_page-1)*nrow_to_print+1
    tmp2=ifelse(print_page*nrow_to_print>nrow(phmdata1),nrow(phmdata1),print_page*nrow_to_print)
    a=pheatmap(phmdata1[seq(tmp1,tmp2),],cluster_rows = F,cluster_cols = F,gaps_col = c(length(significant_cohort_factor)),na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=c(colorRampPalette(corl[1:2])(40),colorRampPalette(corl[2:3])(60)),display_numbers = Cs.mat[seq(tmp1,tmp2),])
  pdf(paste0(gene,".",geneRTtype,".geneCor.pheatmap_pages_",print_page,".pdf"),width = 10,height = 10)
          print(a)
          dev.off()
  }
  }
  setwd(ParentWD)
  cp_report(file.path("28.HallmarkGSVA", outdir))
  message("[step_28_hallmark_gsva] done")
}

step_23_gene_cor <- function() {
  outdir <- paste0("x12.", gene, ".Corr.Heatmap")
  message("[step_23_gene_cor] start")
  setwd(file.path(ParentWD, "23.geneCor"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(limma)

  # for_sale: only immunosuppression-related genes are analyzed here; other pathways are covered by the GSVA step.
  #geneRTtype="immunerepressor"  # gene list type read into geneRT
  #
  #geneRT=read.table(paste0(ParentWD,"/23.geneCor/GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt"),sep="\t",header=F)#emt
  # set geneRT to the gene list to use and run
  #geneRT=read.table(paste0(ParentWD,"/23.geneCor/gene.txt"),sep="\t",header=F)
  #geneRT=read.table(paste0(ParentWD,"/23.geneCor/lipidmetabolism.txt"),sep="\t",header=F)#emt
  #geneRT=read.table(paste0(ParentWD,"/23.geneCor/REACTOME_UNFOLDED_PROTEIN_RESPONSE.txt"),sep="\t",header=F)#UPR

  #geneRTtype="immune"
  #geneRT=read.table(paste0(ParentWD,"/23.geneCor/gene.txt"),sep="\t",header=F)# immune-related genes
  geneRTtype="ICP (Receptor on T&NK cell)"
  geneRT=read.table(paste0(ParentWD,"/23.geneCor/ICPreceptor.txt"),sep="\t",header=F)# immune-related genes

  files=dir(paste0(ParentWD,"/data"))
  files=grep("^symbol.",files,value=T)


  outTab=data.frame()
  corTab=data.frame()
  sameGenes=c()

  pb=txtProgressBar(style=3)# progress bar
  start_time=Sys.time()
  for(i in 1:length(files)){

    CancerType=gsub("^symbol\\.|\\.tsv$","",files[i])
    if( stringr::str_remove(CancerType,"\\.htseq_fpkm") %in% significant_cohort){
    message(paste0("\n","Analysising cancer type: ",CancerType,"\n"))
    setTxtProgressBar(pb,i/length(files))
    rt=readRDS(paste0(ParentWD,"/data/",files[i],".RDS"))

  	rt=as.matrix(rt)
  	rownames(rt)=rt[,1]
  	exp=rt[,2:ncol(rt)]
  	dimnames=list(rownames(exp),colnames(exp))
  	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  	data=avereps(data)

  	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
  	group=sapply(strsplit(group,""),"[",1)
  	group=gsub("2","1",group)
  	data=data[,group==0]

  	sameGenes=intersect(as.vector(geneRT[,1]),row.names(data))
  	data=data[c(sameGenes,gene),]

  	x=as.numeric(data[gene,])
      outVector=data.frame(CancerType)
      corVector=data.frame(CancerType)
  	for(j in sameGenes){
  		y=as.numeric(data[j,])
  		if(sd(y)>0.01){
  			corT=cor.test(x,y)
  			cor=corT$estimate
  			pValue=corT$p.value
  			outVector=cbind(outVector,pValue)
  			corVector=cbind(corVector,cor)
  		}
  		else{
  			outVector=cbind(outVector,pValue=1)
  			corVector=cbind(corVector,cor=0)
  		}
  	}
  	outTab=rbind(outTab,outVector)
  	corTab=rbind(corTab,corVector)
  }
  }
  end_time=Sys.time()
  close(pb)
  run_time=end_time-start_time
  # close progress bar
  colnames(outTab)=c("CancerType",sameGenes)
  outTab$CancerType=as.character(outTab$CancerType)
  outTab$CancerType=stringr::str_remove(outTab$CancerType,"\\.htseq_fpkm")

  write.table(outTab,file=paste0(gene,".",geneRTtype,".geneCor.pvalue.txt"),sep="\t",row.names=F,quote=F)

  colnames(corTab)=c("CancerType",sameGenes)
  corTab$CancerType=as.character(corTab$CancerType)
  corTab$CancerType=stringr::str_remove(corTab$CancerType,"\\.htseq_fpkm")
  write.table(corTab,file=paste0(gene,".",geneRTtype,".geneCor.cor.txt"),sep="\t",row.names=F,quote=F)





  library(reshape2)
  library(RColorBrewer)

  options(stringsAsFactors = F)

  up <- read.table(paste0(gene,".",geneRTtype,".geneCor.pvalue.txt"),sep = "\t",check.names = F,header = T,row.names=1)  
  dn <- read.table(paste0(gene,".",geneRTtype,".geneCor.cor.txt"),sep = "\t",check.names = F,header = T,row.names=1)     
  dn=t(dn)
  up=t(up)

  colVector=c("#AB221F","#3878C1","#FFFADD")

  gene.level <- as.character(rownames(dn)) 
  cancer.level <- as.character(colnames(dn))


  dn.long <- setNames(reshape2::melt(dn), c('Gene', 'Cancer', 'Frequency'))
  dn.long$Categrory <- "DN"
  up.long <- setNames(reshape2::melt(up), c('Gene', 'Cancer', 'Frequency'))
  up.long$Categrory <- "UP"


  dn.long$range <- cut(dn.long$Frequency, 
                       breaks = seq(floor(min(dn.long$Frequency)),
                                    ceiling(max(dn.long$Frequency)),0.01))
  rangeMat1 <- levels(dn.long$range) 
  rbPal1 <- colorRampPalette(colors = c("#2828FE","white","#FF5959"))
  #rbPal1 <- colorRampPalette(colors = c("#0E559D","white","#A93B09"))
  col.vec1 <- rbPal1(length(rangeMat1)); names(col.vec1) <- rangeMat1
  dn.long$color <- col.vec1[as.character(dn.long$range)]


  up.long$range <- cut(up.long$Frequency, breaks = seq(floor(min(up.long$Frequency)),ceiling(max(up.long$Frequency)),0.01)) 
  rangeMat2 <- levels(up.long$range)
  rbPal2 <- colorRampPalette(colors = c("white","grey"))
  col.vec2 <- rbPal2(length(rangeMat2)); names(col.vec2) <- rangeMat2
  up.long$color <- col.vec2[as.character(up.long$range)]


  heatmat <- rbind.data.frame(dn.long,up.long) 
  pdf(paste0(gene,".",geneRTtype,".geneCor.heatmap.pdf"),width = 7,height = 6)
  layout(mat=matrix(c(1,0,1,2,1,0,1,3,1,0),5,2,byrow=T),widths=c(length(cancer.level),2))


  par(bty="n", mgp = c(2,0.5,0), mar = c(5.1, 5.5, 3, 3),tcl=-.25,xpd = T)
  x=as.numeric(factor(heatmat$Cancer,levels = cancer.level))
  y=as.numeric(factor(heatmat$Gene,levels = gene.level))

  plot(1,xlim=c(1,length(unique(x))+1),ylim=c(1,length(unique(y))+1),
       xaxs="i", yaxs="i",xaxt="n",yaxt="n",
       type="n",bty="n",xlab="",ylab="",
       main = "Coexpression across cancer types",cex.main=2)

  for(i in 1:nrow(heatmat)) {
    if(heatmat$Categrory[i]=="DN") polygon(x[i]+c(0,1,1),y[i]+c(0,0,1),col=heatmat$color[i]) 
    if(heatmat$Categrory[i]=="UP") {
      polygon(x[i]+c(0,1,0),y[i]+c(0,1,1),col=heatmat$color[i]) 
      if(heatmat$Frequency[i]<0.001){
        text(x[i]+0.5,y[i]+0.8,"***",cex=0.8)
      }else if(heatmat$Frequency[i]<0.01){
        text(x[i]+0.5,y[i]+0.8,"**",cex=0.8)
      }else if(heatmat$Frequency[i]<0.05){
        text(x[i]+0.5,y[i]+0.8,"*",cex=0.8)
      }
    }
  }

  axis(1,at = sort(unique(x)) + 0.5,labels = cancer.level,lty = 0,las = 2)  
  axis(2,at = sort(unique(y)) + 0.5,labels = gene.level,lty = 0,las = 1)   
  mtext("Cancer types",side = 1,line = 3.5,cex=1.2)  


  par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
  barplot(rep(1,length(col.vec2)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec2)),horiz=TRUE,
          axes = F, col=col.vec2)  # Loss
  axis(4,at=c(1,ceiling(length(col.vec2)/2),length(col.vec2)),c(round(min(up),1),'Pvalue',round(max(up),1)),tick=FALSE)
  par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
  barplot(rep(1,length(col.vec1)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec1)),horiz=TRUE,
          axes = F, col=col.vec1)  # Gain
  axis(4,at=c(1,ceiling(length(col.vec1)/2),length(col.vec1)),c(round(min(dn),1),'Cor',round(max(dn),1)),tick=FALSE)
  dev.off()






  hmpvalue <- read.table(paste0(gene,".",geneRTtype,".geneCor.pvalue.txt"),sep = "\t",check.names = F,header = T,row.names=1)  
  hmdata <- read.table(paste0(gene,".",geneRTtype,".geneCor.cor.txt"),sep = "\t",check.names = F,header = T,row.names=1)     

  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        hmdata[i,j]=0
      }
    }
  }
  hmdata=zz(hmdata)
  hmpvalue=zz(hmpvalue)


  Cs.mat=data.frame(matrix(nrow = nrow(hmpvalue),ncol=ncol(hmpvalue)))
  rownames(Cs.mat)=row.names(hmpvalue)
  colnames(Cs.mat)=colnames(hmpvalue)
  for(i in 1:nrow(hmpvalue)){
    for(j in 1:ncol(hmpvalue)){
      if(hmpvalue[i,j]>=0.05){
        Cs.mat[i,j]='X'
      }else{Cs.mat[i,j]=''}
    }
  }

  corl=c("#244B9E", "white","firebrick")
  a=pheatmap(hmdata,cluster_rows = T,cluster_cols = F,na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
  pdf(paste0(gene,".",geneRTtype,".geneCor.pheatmap_rowcluster.pdf"),width = 12,height = 15)
  print(a)
  dev.off()

  b=pheatmap(hmdata,cluster_rows = F,cluster_cols = F,na_col = "white",border_color = "black",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
  pdf(paste0(gene,".",geneRTtype,".geneCor.pheatmap_nonrowcluster.pdf"),width = 12,height = 15)
  print(b)
  dev.off()

  if(1==2){
    message("running selection pheatmap")
    annotation_row = read.table(paste0(ParentWD,"/23.geneCor/DNAmethylationcluster.txt"),sep="\t",row.names = 1,header = T)
    ann_colors=list(class=c(Methyltransferases='#88CEEA',`Methyl CpG binding proteins`='#D970D5',`Methylcytosine dioxygenases`='#01C689',`Glycosylases`='#F3A462')) # colors for each class
    b=pheatmap(hmdata,cluster_rows = F,annotation_row = annotation_row,annotation_colors=ann_colors,cluster_cols = F,na_col = "white",border_color = "#C7D5D9",cellwidth = 10,cellheight = 10,scale = "none",color=colorRampPalette(corl)(100),display_numbers = Cs.mat)
    pdf(paste0(gene,".",geneRTtype,".geneCor.pheatmap_nonrowcluster.pdf"),width = 12,height = 15)
    print(b)
    dev.off()
  }
  setwd(ParentWD)
  cp_report(file.path("23.geneCor", outdir))
  message("[step_23_gene_cor] done")
}

step_25_gsea <- function() {
  outdir <- paste0("x13.", gene, ".GSEA")
  message("[step_25_gsea] start")
  setwd(file.path(ParentWD, "25.GSEA"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)

  library(limma)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(xlsx)
  library(ggsci)
  library(fgsea)

  gmtGO=read.gmt(paste0(ParentWD,"/25.GSEA/c5.all.v7.1.symbols.gmt"))
  gmtKEGG=read.gmt(paste0(ParentWD,"/25.GSEA/c2.cp.kegg.v7.1.symbols.gmt"))
  gmtGO1=gmtGO
  gmtGO1$term=gsub("GO_","",gmtGO$term)
  gmtGO.list=gmtGO1%>%split(.$term)%>%lapply("[[",2)# for fgsea

  gmtKEGG1=gmtKEGG
  gmtKEGG1$term=gsub("KEGG_","",gmtKEGG$term)
  gmtKEGG.list=gmtKEGG1%>%split(.$term)%>%lapply("[[",2) # for fgsea


  files=dir(paste0(ParentWD,"/data"))
  files=grep("^symbol.*tsv$",files,value=T)
  pb=txtProgressBar(style=3)# progress bar
  start_time=Sys.time()
  for(i in 1:length(files)){

    CancerType=gsub("^symbol\\.|\\.tsv$","",files[i])

  	if( stringr::str_remove(CancerType,"\\.htseq_fpkm") %in% significant_cohort){
  	  message(paste0("Analysising cancer type: ",CancerType,"\n"))
  	  setTxtProgressBar(pb,i/length(files))
  	  rt=readRDS(paste0(ParentWD,"/data/",files[i],".RDS"))


  	rt=as.matrix(rt)
  	rownames(rt)=rt[,1]
  	exp=rt[,2:ncol(rt)]
  	dimnames=list(rownames(exp),colnames(exp))
  	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  	data=avereps(data)

  	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
  	group=sapply(strsplit(group,""),"[",1)
  	group=gsub("2","1",group)
  	data=data[,group==0]

  if(!(gene%in%rownames(data))){next}
  	dataL=data[,(data[gene,]<=median(data[gene,]))]
  	dataH=data[,(data[gene,]>median(data[gene,]))]
  	ncolH=ncol(dataH)
  	ncolL=ncol(dataL)
  	result=c()
  	for(j in 1:nrow(dataH)){result[j]=ifelse(length(grep("\\.",rownames(dataH)[j]))!=0,grep("\\.",rownames(dataH)[j],value = T),NA)}
  	result=result[!is.na(result)]
  	dataH1=dataH[ -which(rownames(dataH) %in%result),]
  	dataL1=dataL[ -which(rownames(dataL) %in%result),]

  	testDF=data.frame(logFC=rep(NA,nrow(dataH1)),pValue=rep(NA,nrow(dataH1)))

  if(all(rownames(dataL1)==rownames(dataH1))){
    for(compare_line in 1:nrow(dataH1)){
      testDF[compare_line,"logFC"]=log2(mean(dataH1[compare_line,],na.rm = T)/(mean(dataL1[compare_line,],na.rm = T)+0.01))
      testDF[compare_line,"pValue"]=wilcox.test(dataH1[compare_line,],dataL1[compare_line,])$p.value
      }
  }     
  	rownames(testDF)=rownames(dataH1)
  	testDF1=testDF[which(testDF$pValue<0.05),]
  	testDF1=testDF1 %>% arrange(desc(logFC))
    log2FC=testDF1$logFC
    names(log2FC)=rownames(testDF1)
  message(paste0("analysising GO gsea in ",CancerType," with gene: ",gene))
      kk=GSEA(log2FC,TERM2GENE=gmtGO, nPerm=1000,pvalueCutoff = 1)
  	kkTab=as.data.frame(kk)
  	kkTab=kkTab[kkTab$pvalue<0.05,]
  	write.xlsx(kkTab,file=paste0(gene,".GSEA.GO.",CancerType,".xlsx"))
  	if((nrow(kkTab)==0)){next}
  	termNum=8
           col_1<-pal_simpsons()(8)# discrete color palette


  	if(nrow(kkTab)>=termNum){
  	printNum=termNum}else{printNum=nrow(kkTab)}

  		gseaplot=gseaplot2(kk, row.names(kkTab)[1:printNum],
  		color=col_1[1:printNum],
  		base_size =14,
  		rel_heights=c(1,0.2,0.4),# relative heights of the panels
  		subplots=1:3,
  		pvalue_table=F,
  		ES_geom="line"
  		)
  		pdf(file=paste0(gene,".GSEA.GO.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".pdf"),width=12,height=10)
  		print(gseaplot)
  		dev.off()

  	# GSEA-style plot	


      fgseaRes <-fgsea(pathways = gmtGO.list,
                stats = log2FC,
                minSize=1,
                maxSize=10000,
                nperm=10000)
      df_fgseaRes = data.frame(lapply(fgseaRes, as.character), stringsAsFactors=FALSE)
      write.table(df_fgseaRes,file = paste0(gene,".fgseaRes.GO.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".txt"),sep="\t",quote = FALSE,row.names = FALSE)  
      topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10),]
      topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), ]
      topPathways <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
      plotG=plotGseaTable(gmtGO.list[topPathways], log2FC, fgseaRes, 
                    gseaParam = 0.5)
      pdf(file=paste0(gene,".fgseaRes.GO.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".pdf"),width=12,height=10)
      print(plotG)
      dev.off()





  	message(paste0("analysising KEGG gsea in ",CancerType," with gene: ",gene))

  	kk1=GSEA(log2FC,TERM2GENE=gmtKEGG, nPerm=1000,pvalueCutoff = 1)
  	kkTab=as.data.frame(kk1)
  	kkTab=kkTab[kkTab$pvalue<0.05,]
  	write.xlsx(kkTab,file=paste0(gene,".GSEA.KEGG.",CancerType,".xlsx"))
  	if((nrow(kkTab)==0)){next}
  	if(nrow(kkTab)>=termNum){printNum=termNum}else{printNum=nrow(kkTab)}
  	gseaplot=gseaplot2(kk1, row.names(kkTab)[1:printNum],
  	color=col_1[1:printNum],
  	base_size =14,
  	rel_heights=c(1,0.2,0.4),# relative heights of the panels
  	subplots=1:3,
  	pvalue_table=F,
  	ES_geom="line")
  	pdf(file=paste0(gene,".GSEA.KEGG.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".pdf"),width=12,height=10)
  	print(gseaplot)
  	dev.off()

  	# GSEA-style plot	


  	fgseaRes <-fgsea(pathways = gmtKEGG.list,
  	                 stats = log2FC,
  	                 minSize=1,
  	                 maxSize=10000,
  	                 nperm=10000)
  	df_fgseaRes = data.frame(lapply(fgseaRes, as.character), stringsAsFactors=FALSE)
  	write.table(df_fgseaRes,file = paste0(gene,".fgseaRes.KEGG.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".txt"),sep="\t",quote = FALSE,row.names = FALSE)  
  	topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10),]
  	topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), ]
  	topPathways <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
  	plotG=plotGseaTable(gmtKEGG.list[topPathways], log2FC, fgseaRes, 
  	                    gseaParam = 0.5)
  	pdf(file=paste0(gene,".fgseaRes.KEGG.",stringr::str_remove(CancerType,"\\.htseq_fpkm"),".pdf"),width=12,height=10)
  	print(plotG)
  	dev.off()
  	}

  }
  end_time=Sys.time()
  close(pb)
  run_time=end_time-start_time
  setwd(ParentWD)
  cp_report(file.path("25.GSEA", outdir))
  message("[step_25_gsea] done")
}

step_29_drugsensitivity <- function() {
  assign("cohort_to_analysis", cohort_focus, envir = globalenv())
  outdir <- paste0("x29.", gene, ".drugsensitivityIC50(pRRophetic)")
  message("[step_29_drugsensitivity] start")
  setwd(file.path(ParentWD, "29.drugsensitivityIC50_pRRophetic"))
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)


  library(limma)
  library(ggpubr)
  library(pRRophetic)
  library(ggplot2)
  set.seed(12345)
  pFilter=0.001 
  allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")# drug list
  Ts=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")# cancer types


  Ts=cohort_to_analysis

  for(cancer in Ts){
  expFile=paste0(ParentWD,"/data/symbol.",cancer,".htseq_fpkm.tsv.RDS")
  rt = readRDS(expFile)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0.5,]
  if(!(gene %in% rownames(data))){
    message(paste0("The gene ",gene," is not found in TCGA ",cancer))
    next
  }
  #colnames(data)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data))
  riskRT=data.frame(data[gene,])
  riskRT=as.matrix(riskRT)
  colnames(riskRT)=gene
  message(paste0("Analysising the correlation of ", gene, " in ",cancer))
  pb=txtProgressBar(style=3)# progress bar
  for(drug in allDrugs){
    message(paste0("Analysising the correlation of ", gene, " in ",cancer," with ", drug, "\n")) 
    setTxtProgressBar(pb,match(drug,allDrugs)/length(allDrugs))
  	senstivity=pRRopheticPredict(data, drug, selection=1,printOutput = F)
  	senstivity=senstivity[senstivity!="NaN"]
  	senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  	sameSample=intersect(row.names(riskRT), names(senstivity))
  	risk=as.data.frame(riskRT[sameSample, gene,drop=F])
  	senstivity=senstivity[sameSample]
  	rt=cbind(risk, senstivity)

    group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
    rt[,gene]=group
  	rt[,gene]=factor(rt[,gene], labels =c("low", "high"),levels=c("low", "high"))
  	type=levels(factor(rt[,gene]))
  	comp=combn(type, 2)
  	my_comparisons=list()
  	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

  	test=wilcox.test(as.formula(paste0("senstivity~",gene)),rt)

  	if(test$p.value<pFilter){

  		boxplot=ggboxplot(rt, x=gene, y="senstivity", fill=gene,
  					      xlab=gene,
  					      ylab=paste0(drug, " Senstivity (Scaled)"),
  					      title=cancer,
  					      palette=c("#244B9E", "firebrick")
  					     )+ 
  			stat_compare_means(comparisons=my_comparisons)+theme(legend.position = "none" ,axis.text.x = element_text(size = 13,colour = "black"),axis.text.y = element_text(size = 11,colour = "black"),axis.title.y = element_text(size = 13,colour = "black"),legend.text = element_text(size = 13,colour = "black"))
  		pdf(file=paste0(gene,".",cancer,".","durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
  		print(boxplot)
  		dev.off()
  	}
  }
  close(pb)
  message(paste0("finishing ",gene," sensitivity analysis in ",cancer))
  }
  setwd(ParentWD)
  cp_report(file.path("29.drugsensitivityIC50_pRRophetic", outdir))
  message("[step_29_drugsensitivity] done")
}

run_imvigor <- function() {
  setwd(file.path(ParentWD, "33.IMvigor210", "IMvigor210"))
  outdir <- paste0("x33.", gene, ".IMvigor210Response")
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
  setwd(outdir)
  #!/usr/bin/env Rscript
  # =====================================================================
  # IMvigor210.geneResponse_for_sale.R - IMvigor210 immunotherapy cohort:
  #   objective response (CR/PR) rates compared between High/Low expression groups
  #
  # Input (parent directory, 33.IMvigor210/IMvigor210/):
  #   ../IMvigor210_cli.txt          clinical table (Best Confirmed Overall Response)
  #   ../IMvigor210_fpkm_symbol.txt   FPKM matrix (genes x 348 samples)
  #
  # Analysis:
  #   1. median split of gene expression (High / Low)
  #   2. response grouping: CR/PR (responder) vs PD/SD (non-responder), NE excluded
  #   3. CR, PR, CR+PR fractions per group; Fisher exact test, High vs Low
  #   4. stacked bar plots: binary (CR/PR vs PD/SD) and four-class (CR/PR/SD/PD)
  #
  # Run mode: set global variable gene, setwd to the output directory,
  #   source("../IMvigor210.geneResponse_for_sale.R")
  # Output written to the working directory:
  #   <gene>.IMvigor210.response.csv
  #   <gene>.IMvigor210.responseBar.{pdf,png}     binary (CR/PR vs PD/SD)
  #   <gene>.IMvigor210.responseBar4.{pdf,png}    four-class (CR/PR/SD/PD)
  # =====================================================================

  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
  })

  if (!exists("gene")) stop("global variable gene (gene symbol) required")
  if (!file.exists("../IMvigor210_cli.txt"))
    stop("input file not found: run from the output directory under 33.IMvigor210/IMvigor210/")

  # ---------------------------- data ----------------------------
  cli <- read.delim("../IMvigor210_cli.txt", sep = "\t", header = TRUE, check.names = FALSE)
  resp_all <- cli[["Best Confirmed Overall Response"]]
  names(resp_all) <- rownames(cli)

  exp_mat <- fread("../IMvigor210_fpkm_symbol.txt", data.table = FALSE)
  rownames(exp_mat) <- exp_mat$genes
  exp_mat <- as.matrix(exp_mat[, -1, drop = FALSE])
  exp_mat <- log2(exp_mat + 1)

  if (!(gene %in% rownames(exp_mat))) stop("gene not present in the expression matrix: ", gene)

  gexpr <- exp_mat[gene, ]
  com <- intersect(names(gexpr), names(resp_all))
  gexpr <- gexpr[com]
  resp <- resp_all[com]

  med <- median(gexpr)
  grp <- factor(ifelse(gexpr > med, "High", "Low"), levels = c("Low", "High"))

  df <- data.frame(sample = com, expr = as.numeric(gexpr), group = grp,
                   Response = resp, stringsAsFactors = FALSE)
  df <- df[df$Response %in% c("CR", "PR", "SD", "PD"), ]
  df$Binary <- ifelse(df$Response %in% c("CR", "PR"), "CR/PR", "PD/SD")

  if (nrow(df) == 0) stop("no samples with a valid response")

  # ---------------------------- proportions and test ----------------------------
  n_by_group <- table(df$group)
  ct4 <- table(factor(df$Response, levels = c("CR", "PR", "SD", "PD")), df$group)
  ct2 <- table(df$Binary, df$group)
  pct4 <- round(100 * prop.table(ct4, margin = 2), 2)
  pct2 <- round(100 * prop.table(ct2, margin = 2), 2)

  ft <- suppressWarnings(fisher.test(ct2))
  p_val <- ft$p.value
  OR <- unname(ft$estimate)

  res_tab <- data.frame(
    gene = gene,
    group = colnames(ct4),
    n = as.integer(n_by_group[colnames(ct4)]),
    CR_n = as.integer(ct4["CR", ]), CR_pct = as.numeric(pct4["CR", ]),
    PR_n = as.integer(ct4["PR", ]), PR_pct = as.numeric(pct4["PR", ]),
    SD_n = as.integer(ct4["SD", ]), SD_pct = as.numeric(pct4["SD", ]),
    PD_n = as.integer(ct4["PD", ]), PD_pct = as.numeric(pct4["PD", ]),
    CRPR_pct = as.numeric(pct2["CR/PR", ]),
    SDPD_pct = as.numeric(pct2["PD/SD", ]),
    odds_ratio = round(OR, 3),
    p_fisher = signif(p_val, 3),
    median_expr = round(med, 3)
  )
  write.csv(res_tab, paste0(gene, ".IMvigor210.response.csv"), row.names = FALSE)

  # ---------------------------- plots ----------------------------
  sci_theme <- theme_classic(base_size = 11) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5))

  save_pub <- function(p, path, w = 5, h = 5.5, dpi = 300) {
    ggsave(paste0(path, ".pdf"), p, width = w, height = h, bg = "white")
    ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  }

  plot_stacked <- function(resp_col, cols, path, w, h) {
    d <- as.data.frame(as.table(prop.table(table(resp_col, df$group), margin = 2)))
    colnames(d) <- c("Response", "Group", "Pct")
    d$label <- sprintf("%.1f%%", 100 * d$Pct)
    p <- ggplot(d, aes(x = Group, y = Pct, fill = Response)) +
      geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.3) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5),
                size = 3.2, colour = "white") +
      scale_fill_manual(values = cols) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                         labels = scales::percent_format()) +
      labs(x = paste0(gene, " expression"), y = "Proportion",
           fill = "Response", title = gene) +
      annotate("text", x = 1.5, y = 1.02,
               label = paste0("Fisher exact P = ", signif(p_val, 3)),
               size = 3.2) +
      coord_cartesian(ylim = c(0, 1.08)) +
      sci_theme
    save_pub(p, path, w, h)
  }

  plot_stacked(df$Binary, c("CR/PR" = "#D6604D", "PD/SD" = "#4393C3"),
               paste0(gene, ".IMvigor210.responseBar"), 5, 5.5)
  plot_stacked(factor(df$Response, levels = c("CR", "PR", "SD", "PD")),
               c("CR" = "#B2182B", "PR" = "#EF8A62", "SD" = "#92C5DE", "PD" = "#2166AC"),
               paste0(gene, ".IMvigor210.responseBar4"), 5.5, 5.5)

  message("IMvigor210 response analysis done: ", gene,
          " (n High = ", res_tab$n[2], ", Low = ", res_tab$n[1],
          ", CR+PR High = ", res_tab$CRPR_pct[2], "% vs Low = ", res_tab$CRPR_pct[1],
          "%, P = ", signif(p_val, 3), ")")
  setwd(ParentWD)
  cp_report(file.path("33.IMvigor210", "IMvigor210", outdir))
}


# =====================================================================
# 4. Main
# =====================================================================
.rd <- file.path(ParentWD, "report", gene)
if (!dir.exists(.rd)) dir.create(.rd, recursive = TRUE, showWarnings = FALSE)

.steps <- list(
  step_06_diff, step_08_survival,
  step_09_cox_os, step_10_cox_dss, step_11_cox_dfi, step_12_cox_pfi,
  step_16_tmb_radar, step_18_msi_radar,
  step_20_estimate, step_21_cibersort, step_22_cibersort_cor,
  step_26_immunefunction, step_27_timer2, step_28_hallmark_gsva,
  step_23_gene_cor, step_25_gsea, step_29_drugsensitivity
)

for (.st in .steps) {
  tryCatch(.st(),
           error = function(e) {
             message("[FAIL] ", conditionMessage(e))
             try(setwd(ParentWD), silent = TRUE)
           })
  if (identical(.st, step_08_survival)) {
    odir <- file.path(ParentWD, "08.survival", paste0("2.", gene, ".survivalplot"))
    sig <- c()
    for (f in grep("\\.pdf$", list.files(odir), value = TRUE)) {
      tmp <- unlist(strsplit(f, "\\."))
      sig <- union(sig, tmp[length(tmp) - 1])
    }
    assign("significant_cohort", sig, envir = globalenv())
    message("[main] significant cohorts: ", paste(sig, collapse = ", "))
  }
}

if (dir.exists(file.path(ParentWD, "33.IMvigor210", "IMvigor210"))) {
  tryCatch(run_imvigor(),
           error = function(e) { message("[FAIL] IMvigor: ", conditionMessage(e)); try(setwd(ParentWD), silent = TRUE) })
}

message("======== pipeline finished, results in report/", gene, "/ ========")
