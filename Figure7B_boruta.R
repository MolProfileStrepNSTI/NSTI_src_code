##################################################
#                                                #
#    Boruta feature selection for Figure 7B      #
#                                                #
##################################################

# Autoinstall missing packages
# src: https://gist.github.com/benmarwick/5054846#gistcomment-2613621
list.of.packages = c("Boruta") # replace xx and yy with package names
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)


# data import
data_file <- "data/fig7/train_data_all_markers.tsv"
data <- read.csv(data_file, sep='\t')

columns <- c("sample", "outcome", "IL_1b", "IL_6", "TNF_a", "CXCL9",
             "CXCL10", "CXCL11", "IFN_lambda", "IL_8", "IL_12p70", "IFN_a", 
             "IL_28", "GM_CSF", "IFN_b", "IL_10", "IFN_g")
colnames(data) <- columns

# data cleaning
data$outcome <- toupper(data$outcome)
# remove samples with missing values
data <- na.omit(data)
# exclude Healthy Patient Samples from Dataset
data <- data[data$outcome != 'HEALTHY', ]

rownames(data) <- data$sample
data = subset(data, select=-sample)
data$outcome <- factor(data$outcome)

data$outcome <- as.integer(data$outcome == 'STREPTOCOCCUS PYOGENES')

# Boruta feature selection
boruta.chemokine <- Boruta(outcome~., data=data, doTrace=2, ntree=100)

resultDir <- 'results/figure7'
dir.create(resultDir, recursive=TRUE)

png("results/figure7/Fig7B_Boruta.png", width = 640, height = 480)

par(mar=c(8, 4.1, 4.1, 2.1))
plot(boruta.chemokine, las=2, xlab='', main='Boruta - Vanilla')
dev.off()

print(boruta.chemokine)

# apply RoughFix if needed
if('Tentative' %in% boruta.chemokine$finalDecision) {
  print('---- Fixing Tentative Features ----')
  borTenFix.chemokine = TentativeRoughFix(boruta.chemokine)
  
  png("results/figure7/Fig7B_Boruta_tentative_fixed.png", width = 640, height = 480)
  par(mar=c(8, 4.1, 4.1, 2.1))
  plot(borTenFix.chemokine, las=2, xlab='', main='Boruta - Tentative Fixed')
  dev.off()
  
  print(borTenFix.chemokine)
}