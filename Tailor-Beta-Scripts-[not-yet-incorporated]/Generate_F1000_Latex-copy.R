## ----loadLibs, eval=FALSE--------------------------------------------------
  source("http://www.bioconductor.org/biocLite.R")
  biocLite('BiocWorkflowTools')

## ----RStudioNew, echo = FALSE, out.width='90%', fig.align='center', fig.cap="Creation of a new article is integrated into RStudio.  The F1000Research template can be accessed via the 'new R Markdown document' menu option."----
knitr::include_graphics("Rstudio_newdoc.png")

## ----createDraft-----------------------------------------------------------
tmp_dir <- tempdir()
setwd(tmp_dir)

rmarkdown::draft(file = "MyArticle.Rmd",
                 template = "f1000_article",
                 package = "BiocWorkflowTools",
                 edit = FALSE)

## ----newFiles, echo = FALSE, out.width='90%', fig.align='center', fig.cap="The complete set of files.  After 'knitting' both the LaTeX source file and PDF documents can be found alongside your Rmarkdown file."----
knitr::include_graphics("new_files.png")

## ----render, eval = FALSE--------------------------------------------------
#  rmd_file <- file.path(tmp_dir, 'MyArticle', 'MyArticle.Rmd')
#  rmarkdown::render(input = rmd_file)

## ----upload, eval=FALSE----------------------------------------------------
  library(BiocWorkflowTools)
#  workflow_dir <- file.path(tmp_dir, 'MyArticle')
#  uploadToOverleaf(files = workflow_dir)

## ----sessionInfo, eval=TRUE------------------------------------------------
sessionInfo()

## ---- include_logo, out.width="50%", echo=FALSE----------------------------
knitr::include_graphics("denbi_logo.jpg")

