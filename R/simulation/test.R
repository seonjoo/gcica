.libPaths("/ifs/scratch/msph/LeeLab/softwares/R/hpc")
library(ggplot2)

dat=data.frame(x=rnorm(100))
png('example.png')
ggplot(dat,aes(x))+geom_histogram()
dev.off()
save.image('test.Rdata')