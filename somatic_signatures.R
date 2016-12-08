library(SomaticSignatures)
vcf <- readVcf("HPV16_tumor_snv_files/CESC-US-16_withvirus_chr_merge.vcf.gz", genome="hg19")
vcf2 <- readVcf("HPV16_tumor_snv_files_novirus/CESC-US-16_novirus_chr_merge.vcf.gz", genome="hg19")
vr <- as(vcf, "VRanges")
vr2 <- as(vcf2, "VRanges")


ctx = mutationContext(vr, BSgenome.Hsapiens.UCSC.hg19)
sampleNames(ctx) = "sample1"
ctx2 = mutationContext(vr2, BSgenome.Hsapiens.UCSC.hg19)
sampleNames(ctx2) = "sample2"


m = motifMatrix(ctx, group = "sampleNames")
m2 = motifMatrix(ctx2, group = "sampleNames")

mut_comb <- cbind(m,m2)

sig <- identifySignatures(comb, 2)
