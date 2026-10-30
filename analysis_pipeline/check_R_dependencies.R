pkgs<-c('mice','sandwich','metafor','numDeriv')
for(p in pkgs){if(!requireNamespace(p,quietly=TRUE))stop(paste('Install required R package:',p));cat(p,as.character(packageVersion(p)),'\n')}
