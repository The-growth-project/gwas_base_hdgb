args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

df <- read.csv(input_file, sep = "", header = TRUE)
head(df)

if (!all(c("ID", "LOG10P") %in% colnames(df))) {
    stop("The columns ID and LOG10P are not present")
}

df_out <- data.frame(SNP = df$ID, CHR = df$CHROM, POS = df$GENPOS , A1 = df$ALLELE1,A2 = df$ALLELE0, MAF = df$A1FREQ,  BETA = df$BETA, SE = df$SE,  P = 10^(-df$LOG10P))
df_out <- df_out[df_out$MAF>=0.001,]
df_out <- df_out[df_out$P<5e-08,]
head(df_out)

write.table(df_out, file = output_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
