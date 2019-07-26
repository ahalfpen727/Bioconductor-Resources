## ---- echo=FALSE, comment=NA---------------------------------------------
library(httr)
cat(x <- content(GET("https://raw.githubusercontent.com/dtenenba/bioc_docker/master/tracker.bioconductor.org/fig.yml")))

