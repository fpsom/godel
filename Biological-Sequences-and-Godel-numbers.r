library(seqinr)
library(MASS)
library(ggplot2)
library(ggpubr)

sieve <- function(n) {
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
  {
    primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
    last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
  }
  which(primes)
}

seqList <- function(x, type) {
  if (type == "DNA") {
    as.list(read.fasta(x,
                       as.string = FALSE,
                       seqtype="DNA",
                       seqonly = FALSE,
                       strip.desc = TRUE))
  }
  else {
    if (type == "AA") {
      as.list(read.fasta(x,
                         as.string = FALSE,
                         seqtype="AA",
                         seqonly = FALSE,
                         strip.desc = TRUE))
    }
    else {
      NULL
    }
  }
}

createRandomSequencesBasedOnDistr <- function(count, length, prob=c(0.25,0.25,0.25,0.25), fileNameRandSeqs) {
  
  sink(fileNameRandSeqs)
  set.seed(seedValuesList[1]) 
  for (i in 1:count) {
    cat(">Seq", i, "\n", sep = "")
    seqX <- sample(c("A","C","G","T"), length, rep=TRUE, prob)
    cat(paste(seqX,collapse=""), "\n", sep = "")
  }
  sink()
}

createRandomSequenceValues <- function(seedList, type) {
  
  if ( type == "DNA" ) {
    results <- matrix(, nrow = length(seedList), ncol = 4)
  }
  else if ( type == "AA" ) {
    results <- matrix(, nrow = length(seedList), ncol = 23)
  }
  
  for (i in 1:length(seedList)) {
    set.seed(seedValuesList[i])
    if ( type == "DNA" ) {
      results[i,] <- sample(seq(from = 1, to = 4, by = 1), size = 4, replace = FALSE)
    }
    else if ( type == "AA" ) {
      results[i,] <- sample(seq(from = 1, to = 23, by = 1), size = 23, replace = FALSE)
    }
  }
  return(results)
}

assignSets <- function(randSequenceValues, type) {
  
  stringValues <- ""
  
  if (type == "DNA") {
    #DNA Set
    cat("Sequence Values: (")
    
    for (i in 1:4) {
      stringValues <- paste(stringValues, randSequenceValues[i], ", ")
    }
    cat(stringValues)
    cat(")\n")
    assign("a", randSequenceValues[1], envir = .GlobalEnv)
    assign("t", randSequenceValues[2], envir = .GlobalEnv)
    assign("g", randSequenceValues[3], envir = .GlobalEnv)
    assign("c", randSequenceValues[4], envir = .GlobalEnv)
  }
  else if ( type == "AA") {
    # AA Set 1
    cat("Sequence Values: (")
    for (i in 1:23) {
      stringValues <- paste(stringValues, randSequenceValues[i], ", ")
    }
    cat(stringValues)
    cat(")\n")
    assign("A", randSequenceValues[1], envir = .GlobalEnv) # alanine, ala
    assign("R", randSequenceValues[2], envir = .GlobalEnv) # arginine, arg
    assign("N", randSequenceValues[3], envir = .GlobalEnv) # asparagine, asn
    assign("D", randSequenceValues[4], envir = .GlobalEnv) # aspartic acid, asp
    assign("B", randSequenceValues[5], envir = .GlobalEnv) # sparagine or aspartic acid, asx
    assign("C", randSequenceValues[6], envir = .GlobalEnv) # cysteine, cys
    assign("E", randSequenceValues[7], envir = .GlobalEnv) # glutamic acid, glu
    assign("Q", randSequenceValues[8], envir = .GlobalEnv) # glutamine, gln
    assign("Z", randSequenceValues[9], envir = .GlobalEnv) # glutamine or glutamic acid, glx
    assign("G", randSequenceValues[10], envir = .GlobalEnv) # glycine, gly
    assign("H", randSequenceValues[11], envir = .GlobalEnv) # histidine, his
    assign("I", randSequenceValues[12], envir = .GlobalEnv) # isoleucine, ile
    assign("L", randSequenceValues[13], envir = .GlobalEnv) # leucine, leu
    assign("K", randSequenceValues[14], envir = .GlobalEnv) # lysine, lys
    assign("M", randSequenceValues[15], envir = .GlobalEnv) # methionine, met
    assign("F", randSequenceValues[16], envir = .GlobalEnv) # phenylalanine, phe
    assign("P", randSequenceValues[17], envir = .GlobalEnv) # proline, pro
    assign("S", randSequenceValues[18], envir = .GlobalEnv) # serine, ser
    assign("T", randSequenceValues[19], envir = .GlobalEnv) # threonine, thr
    assign("W", randSequenceValues[20], envir = .GlobalEnv) # tryptophan, trp
    assign("Y", randSequenceValues[21], envir = .GlobalEnv) # tyrosine, tyr
    assign("V", randSequenceValues[22], envir = .GlobalEnv) # valine, val
    assign("X", randSequenceValues[23], envir = .GlobalEnv) # undetermined
  }
  
  stringValues
}

godelStatistics <- function(x) {
    
  if (logOutput) {
    cat(paste("  Min  : ", summary(x)[1], "\n"))
    cat(paste("1st Qu.: ", summary(x)[2], "\n"))
    cat(paste("Median : ", summary(x)[3], "\n"))
    cat(paste(" Mean  : ", summary(x)[4], "\n"))
    cat(paste("3rd Qu.: ", summary(x)[5], "\n"))
    cat(paste("  Max  : ", summary(x)[6], "\n"))
    cat(paste("St. Dev: ", sd(x), "\n"))
  }
  
  results <- list(minG = summary(x)[1], firstQG = summary(x)[2],
       medianG = summary(x)[3], meanG = summary(x)[4], 
       thirdQG = summary(x)[5], maxG = summary(x)[6],
       stdG = sd(x))
  
  return(results)
}

primes <- sieve(20000) # length of primes should be >= max sequence length
logOutput <- FALSE     # debug output
numberOfPoints <- 4;   # number of different assignments of letters for Godel numbers
replicate <- FALSE     # if we need to set specific seed numbers

if (replicate == FALSE) {
  seedValuesList <- sample(seq(from = 1, to = 1000, by = 1), size = numberOfPoints, replace = FALSE)
  for (i in 1:numberOfPoints) {
    cat(paste(seedValuesList[i], ", "))
  }
} else {
  seedValuesList <- as.numeric(unlist(read.csv(file="seeds.csv", header = FALSE)))
  numberOfPoints <- length(seedValuesList)
}

type <- "DNA"
# The following values should execute within 6 minutes with a good enough resolution
seqLengthLimit <- 361
numberOfArtificialSeqs <- 821     # "A","C","G","T"
# amend the following line for your particular distribution on nucloetide presence.
createRandomSequencesBasedOnDistr(numberOfArtificialSeqs, seqLengthLimit, c(0.25,0.25,0.25,0.25), "data/artificialSeqs.fasta")
ompGene.list <- seqList("data/artificialSeqs.fasta", type)
# ompGene.list <- seqList("data/realSeqs.fasta", type)
ompGene.list.raw <- ompGene.list

randValues <- createRandomSequenceValues(seedValuesList, type)

randValues

# Do not run this cell if you want to check for all nucleotides and not just for A
# randValues[1,] <- c(3, 2, 1, 4)
randValues

ompGene.list <- ompGene.list.raw[which(getLength(ompGene.list.raw) >= 1)]
sizeExp <- length(ompGene.list)
selectedSequences <- c(1:sizeExp)

godelValuePoints <- data.frame(matrix(0, ncol = numberOfPoints, nrow = sizeExp))
namesList <- list()
for (i in 1:numberOfPoints) {
  namesList[i] <- paste('godel_log_pos', i, sep = '')
}
godelValuePoints <- setNames(godelValuePoints, namesList)
godelValuePoints$seqNames <- unlist(attributes(ompGene.list)$name)

stringAssignValues <- vector(mode="character", length=numberOfPoints)

for (indexPos in 1:numberOfPoints) {
  stringAssignValues[indexPos] <- assignSets(randValues[indexPos,], type)
  
  godel.value.exp <- list()
  godel.value.log <- list()
  
  for (indexSeq in selectedSequences) {
    
    godel.value.exp[[indexSeq]] <- 1
    godel.value.log[[indexSeq]] <- 0
    
    if (logOutput) {
      cat("Sequence length  ")
      cat(indexSeq)
      cat("  :")
      cat(length(ompGene.list[indexSeq][[1]]))
      cat("\n")
    }

    for (i in 1:length(ompGene.list[indexSeq][[1]])) {
      
      prime <- as.numeric(primes[i])
      alpha <- as.numeric(get(as.character(ompGene.list[[indexSeq]])[i]))
      
      godel.value.exp[[indexSeq]] <- godel.value.exp[[indexSeq]] * prime ** alpha
      godel.value.log[[indexSeq]] <- godel.value.log[[indexSeq]] + alpha*log(prime)
    }
    
  }
  
  godelValuePoints[, paste('godel_log_pos', indexPos, sep = '')] <- unlist(godel.value.log)
}

statsPos <- data.frame(matrix(0, ncol = 7, nrow = numberOfPoints))
statsPos <- setNames(statsPos, c("minG", "firstQG", "medianG", "meanG", "thirdQG", "maxG", "stdG"))
for (indexPos in 1:numberOfPoints) {
  statsPos[indexPos, ] <- godelStatistics(godelValuePoints[, paste('godel_log_pos', indexPos, sep = '')])
}

statsPos

P1 <- sum(log(primes[1:seqLengthLimit]))
P2 <- sum((log(primes[1:seqLengthLimit]))^2)

theoreticalMeanEqual <- P1*2.5
theoreticalStdEqual <- sqrt(P2*1.25)

theoreticalMeanEqual
theoreticalStdEqual


para <- matrix(nrow=numberOfPoints, ncol = 2)
for (indexPos in 1:numberOfPoints) {
  fit <- fitdistr(godelValuePoints[, indexPos], "normal")
  para[indexPos, 1] <- fit$estimate["mean"]
  para[indexPos, 2] <- fit$estimate["sd"]
}


statsPos$stringAssignValues <- as.factor(stringAssignValues)
statsPos$fit_estimate <- ((para[,1]-theoreticalMeanEqual)/para[,1])*100
statsPos$fit_sd <- ((para[,2]-theoreticalStdEqual)/para[,2])*100



mean(statsPos$fit_estimate)
mean(statsPos$fit_sd)


indexPos = 1

binwidthPlot = 5

ggplot(godelValuePoints, aes(x=godelValuePoints[, indexPos]), environment = environment()) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=binwidthPlot)+
  geom_density(alpha=.2, fill="#FF6666") +
  stat_function(fun = dnorm, args = list(mean = para[indexPos, 1], sd = para[indexPos, 2]), color = "darkred", size = 2, linetype = "dotdash") +
  labs(title= "Histogram and Density plot of Godel numbers",x="Godel numbers", y = "Density")+
  scale_color_brewer(palette="Accent") + 
  theme_minimal()



p1 <- ggplot(statsPos, aes(x = reorder(stringAssignValues, fit_estimate), y = fit_estimate, group=1)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Difference between theoretical and observed value of Mean (%)",x="Assignments", y = "Difference in mean (%)")

p2 <- ggplot(statsPos, aes(x = reorder(stringAssignValues, fit_estimate), y = fit_sd, group=1)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Difference between theoretical and observed value of standard deviation (%)",x="Assignments", y = "Difference in  standard deviation (%)")

ggarrange(p1 + rremove("x.text"), p2,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)


# a t g c
lapply(which(statsPos$maxG == min(statsPos$maxG)),
       function(x) assignSets(randValues[x,], "DNA"))
# artificial dataset: 4 , 1 , 3 , 2 (min of meanG)
# artificial dataset: 4 , 2 , 1 , 3 (min of minG)
# artificial dataset: 2 , 1 , 3 , 4 (min of maxG)

# a t g c
lapply(which(statsPos$maxG == max(statsPos$maxG)),
       function(x) assignSets(randValues[x,], "DNA"))
# artificial dataset: 1 , 4 , 2 , 3 (max of meanG)
# artificial dataset: 4 , 3 , 2 , 1 (max of minG)
# artificial dataset: 1 , 3 , 4 , 2 (max of maxG)


