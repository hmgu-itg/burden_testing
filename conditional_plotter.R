#!/software/R-3.3.0/bin/R CMD BATCH --vanilla --no-save

library("grDevices")

# Assign randon colors:
colorize = function (x){
    if(is.na(x)){
        col = "black"
    }
    else {
        col = palette(rainbow(max_color))[x]
    }
    if( length(col) == 0 ){col = "black"}    
    return(adjustcolor(col, alpha.f = 0.4))
}

# fileName="TMEM40_conditionals.tsv_ext"
args = commandArgs(TRUE)
fileName = args[1]

cat(getwd())

# Reading dataframe:
df = read.table(fileName, sep="\t", header=TRUE, as.is=TRUE)

# Calculate log10 conditional p-values:
df$logPval = log10(as.numeric(df$Cond_pVal))

# Calculate genomic coordinates:
df$coordinates = as.vector(sapply(df$Variant, function(x) as.numeric(strsplit(x,":")[[1]][2]), simplify=TRUE))

# Assingning colors to LD blocks:
max_color = max(df$LD, na.rm=T)
df$colors = sapply(df$LD, colorize, simplify=T)

# Calculate burden p-values:
bPval = unique(df$Burden_pVal)[1]

# Gene_name:
geneName = strsplit(fileName, "_")[[1]][1]

# Generate labels:
x_label = paste("Genomic coordinates on",strsplit(df$Variant[1], ":")[[1]][1])
main_test = paste(geneName, " - conditional and single-point p-values")

# Initialize plot:
pdf(paste(geneName,"_conditionals.pdf", sep=''), 12, 6)
par(mfrow = c(2, 1), oma = c(1,2,0,2))

par(mar=c(0, 4, 4, 1))
plot(y = -log10(df$SpPval),
     x = df$coordinates,
     ylab="-log10(SP p-value)",
     col=df$colors,
     main=main_test,
     ylim=c(-0.2, -log10(bPval)*1.1),
     pch=19,
     xaxt='n')
abline(h=-log10(bPval), col="firebrick")

par(mar=c(4, 4, 0, 1))
plot(y = df$logPval,
     x = df$coordinates,
     ylab="log10(burden p-value)",
     xlab = x_label,
     ylim=c(log10(bPval)*1.1, 0),
     col=df$colors, pch=19)
abline(h=log10(bPval), col="firebrick")
dev.off()
