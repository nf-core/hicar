#!/usr/bin/env Rscript
pkgs <- commandArgs(trailingOnly = TRUE)
print(pkgs)
# set libPath to pwd
lib <- .libPaths()
if(file.access(lib[1], mode=2)!=0){
    pwd <- getwd()
    pwd <- file.path(pwd, "lib")
    dir.create(pwd)
    .libPaths(c(pwd, lib))
}

if(length(pkgs)>0){
    while(!requireNamespace("BiocManager", quietly = TRUE)){
        tryCatch(
            install.packages("BiocManager",
                            repos = "https://cloud.r-project.org/",
                            quiet = TRUE),
            error = function(.e){
                message("retry package installation.")
            }
        )

    }
    pkgs <- pkgs[pkgs %in% BiocManager::available() | grepl("\\/", pkgs)]

    if(any(grepl("\\/", pkgs))){
        pkgs <- c("remotes", pkgs)
    }
    Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
    getPkg <- function(pkgs){
        for(pkg in pkgs){
            retryCount <- 0
            while(!requireNamespace(pkg, quietly = TRUE) && retryCount<3){
                retryCount <- retryCount+1
                tryCatch(
                    if(grepl("\\/", pkg)){
                        remotes::install_github(pkg, upgrade = FALSE, quite = FALSE)
                    }else{
                        BiocManager::install(pkg, update = FALSE, ask = FALSE)
                    },
                    error = function(.e){
                        message("retry package installation.")
                        Sys.sleep(60)
                    }
                )
            }
        }
    }
    getPkg(pkgs)
}
