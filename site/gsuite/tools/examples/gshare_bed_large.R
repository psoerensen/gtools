# Bounded packed-BED execution check using an existing gsim 50K-marker panel.
# Run from gsuite. This tests 600 selected markers, not genome-wide scale.
.libPaths(c("build/r-library",.libPaths()))
library(gsuite)
root <- "build/examples/simulated-genomics"
record <- file.path(root,"data.rds")
if(!file.exists(record)) stop("The shared gsim panel is not prepared: ",record)
x <- readRDS(record)$value
markers <- unlist(x$Glist$rsids,use.names=FALSE)[1:600]
frequency <- setNames(unlist(x$Glist$af,use.names=FALSE)[1:600],markers)
ids <- x$Glist$ids
response <- x$simulation$Y[,1]
stopifnot(length(ids)==length(response),length(markers)==600L,
  all(is.finite(response)),all(is.finite(frequency)&frequency>0&frequency<1))
selected <- list(A=seq_len(500),B=501:2000)
out <- "build/examples/gshare-bed-large"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fit <- function(site,transfer=NULL) {
  i <- selected[[site]]
  Glist <- gprep(bedfiles=file.path(root,"genotypes.bed"),
    ids=ids[i],rsids=markers)
  prior <- if(is.null(transfer))
    list(residual_variance=1,effect_variance=.04,inclusion_probability=.1) else
    list(residual_variance=1)
  gshare(data=data.frame(id=ids[i],trait=response[i]-mean(response[i])),
    Glist=Glist,reference_frequency=frequency,site=site,
    trait="trait",units="centered simulation units",
    covariates="site-centered phenotype; no other covariates",
    prior=prior,transfer=transfer,disjoint=!is.null(transfer),
    control=list(burnin=100L,sampling_sweeps=300L,seeds=31L,threads=1L))
}
source <- fit("A")
message <- gshare_transfer(source,.7,.0204,probability_floor=1e-6)
write_gshare(message,file.path(out,"source.gshare"),overwrite=TRUE)
receiver <- fit("B",read_gshare(file.path(out,"source.gshare")))
peak_command <- paste0("(Get-Process -Id ",Sys.getpid(),").PeakWorkingSet64")
peak_bytes <- as.numeric(system2("powershell.exe",
  c("-NoProfile","-Command",shQuote(peak_command)),stdout=TRUE))
summary <- data.frame(markers=length(markers),source_n=length(selected$A),
  receiver_n=length(selected$B),source_seconds=source$seconds,
  receiver_seconds=receiver$seconds,process_peak_mib=peak_bytes/1024^2,
  source_zero_sampled=sum(source$effects$pip==0),
  source_positive_conditional=sum(source$effects$conditional_inclusion>0))
write.csv(summary,file.path(out,"execution.csv"),row.names=FALSE)
print(summary)
