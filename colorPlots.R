pd <- data.frame(PC1=H[,1], PC2=H[,2], PC3=H[,3],
                 population=as(KG_SAMPLE_POPL_CODE$Population[], "vector"))

popnames <- list(
  ASW='African-Amer',
  CEU='European',
  CHB='Asian',
  CHS='Asian',
  CLM='Colombian',
  FIN='European',
  GBR='European',
  IBS='European',
  JPT='Asian',
  MXL='Mexican-Amer',
  PUR='Puerto Ricans',
  TSI='European',
  YRI='African',
  LWK='African'
)
for (i in 1:length(popnames)) {
  levels(pd$population)[levels(pd$population) == names(popnames)[i]] <- popnames[[names(popnames)[i]]]
}
pd$population <- factor(pd$population,levels(pd$population)[c(13,2,5,7,18,15,6)])


cbbPalette <- c("#E69F00", "#CC79A7", "#D55E00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#FF0000")

ggplot(pd, aes(PC1,PC2,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))
