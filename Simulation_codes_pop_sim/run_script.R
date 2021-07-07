rm(list=ls())
library(foreach)
library(doMC)
n <- 100
registerDoMC(4)

seeds <- sample(c(1:100000000), n,  replace=FALSE)
#seeds<-c(1:1000000000, n, replace=FALSE)
#print(typeof(seeds))

print(length(seeds))

for (i in 1:length(seeds)){
    print(seeds[i])
    }


system("make")

curDir <- getwd()

name <- ""

foreach(i =  1:n)%dopar%{
    seed<-seeds[i]
    cmds <- paste(curDir, "/pop_sim -N 100 -gens 100000 -linear -gom_time_dep -deterministic -non_herited -non_branching -seed ", seeds[i],  sep="")
    setwd(file.path(curDir, name))
    dirName <- paste("results", i, sep="")
    workingdir<-file.path(paste(curDir, "/", name, sep=""), dirName)
    dir.create(workingdir, showWarnings=FALSE)

    setwd(workingdir)
    ##print(c(workingdir, "*****"))
    #system(paste("cp ", curDir, "/spat*.txt .", sep=""))

    cmd <- cmds
    print(cmd)
    system(cmd)
    getwd()
 }
