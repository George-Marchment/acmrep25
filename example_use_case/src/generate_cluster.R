
args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1],header = TRUE, sep = ",")
pdf(file = "plot.pdf") 
plot(data$x, data$y,
     col = factor(data$risk))
dev.off()

