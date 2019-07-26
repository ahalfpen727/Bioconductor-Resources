# 1. Before you upgrade, build a temp file with all of your old packages.
tmp <- installed.packages()
installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpkgs, file="installed_old.rda")

# 2. Install the new version of R and let it do it’s thing.
# 3. Once you’ve got the new version up and running, reload the saved packages and re-install them from CRAN.

load("installed_old.rda")
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
missing <- setdiff(installedpkgs, installedpkgs.new)
install.packages(missing)
update.packages()

# Note: If you had any packages from BioConductor, you can update those too!

source("http://bioconductor.org/biocLite.R")
chooseBioCmirror()
biocLite()
load("installed_old.rda")
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
missing <- setdiff(installedpkgs, installedpkgs.new)
for (i in 1:length(missing)) biocLite(missing[i])

R_LIBS_USER<-file.path("~/R/x86_64-pc-linux-gnu-library/3.4")
fs::dir_copy(dirname("~/R/x86_64-pc-linux-gnu-library/3.6/"), new_path="~/R/x86_64-pc-linux-gnu-library/3.4")
fs::dir_create(Sys.getenv("R_LIBS_USER"))
# Create a new directory for the version of R
fs::dir_create("~/Library/R/3.4/library")
update.packages()
usethis::edit_r_environ()
library(installr)
copy.packages.between.libraries(from, to, ask = FALSE, keep_old = TRUE,
                                do_NOT_override_packages_in_new_R = TRUE)

# Re-start R so the .libPaths are updated
chooseCRANmirror()
library(devtools)
devtools::edit_r_environ()
installr::updateR()
# Lookup what packages were in your old package library
pkgs <- fs::path_file(fs::dir_ls("~/R/x86_64-pc-linux-gnu-library/3.5/"))
devtools::update_packages("~/R/x86_64-pc-linux-gnu-library/3.5/")
