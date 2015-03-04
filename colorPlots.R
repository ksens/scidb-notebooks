H <- SCVsvd$u[,0:4][] %*% diag(SCVsvd$d[0:4][])
pd <- data.frame(PC1=H[,1], PC2=H[,2], PC3=H[,3], PC4=H[,4], PC5=H[,5],
                 population=as(KG_SAMPLE_POPL_CODE$Population[], "vector"))

popnames <- list(
  ACB='Caribbean',
  ASW='African-Amer',
  BEB='South Asian',    # Bengali in Bangladesh
  CEU='European',
  CHB='East Asian',
  CDX='East Asian',    # Chinese Dai in Xishuangbanna
  CHS='East Asian',
  CLM='South American',
  ESN='African',        # Esan in Nigeria
  FIN='European',
  GBR='European',
  GIH='South Asian',    # Gujarati Indian in Houston
  GWD='African',        # Gambian in Western Division
  IBS='European',
  ITU='South Asian',    # Indian Telugu in the UK
  JPT='East Asian',
  KHV='East Asian',    # Kinh in Ho Chi Minh City
  LWK='African',
  MSL='African',        # Mende in Sierra Leone
  MXL='Mexican-Amer',
  PEL='South American', # Peruvian in Lima
  PJL='South Asian',    # Punjabi in Lahore, Pakistan
  PUR='Caribbean',      # 'Puerto Ricans'
  STU='South Asian',
  TSI='European',
  YRI='African'
)
for (i in 1:length(popnames)) {
  levels(pd$population)[levels(pd$population) == names(popnames)[i]] <- popnames[[names(popnames)[i]]]
}
#pd$population <- factor(pd$population,levels(pd$population)) #[c(13,2,5,7,18,15,6)])


cbbPalette <- c("#E69F00", "#CC79A7", "#D55E00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#FF0000",
                "maroon")

## Plot PC1 vs PC2
ggplot(pd, aes(PC1,PC2,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))

## Plot PC1 vs PC3
ggplot(pd, aes(PC1,PC3,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))

## Plot PC2 vs PC3
ggplot(pd, aes(PC2,PC3,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))




