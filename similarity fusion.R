#数据来源
############
#相似度融合：在两个disease或者lncRNA相似度矩阵的第[i,j]个位置取最大值
##################
Jaccard_disease<-read.csv(file = "D:/科研论文/实验材料/disease_Jacad.csv",header=F)
Jaccard_lncRNA<-read.csv(file = "D:/科研论文/实验材料/lncRNA_Jacad.csv",header=F)
GIP_lncRNA<-read.csv(file = "D:/科研论文/实验材料/GIP_lncRNA.csv",header=F)
GIP_disease<-read.csv(file = "D:/科研论文/实验材料/GIP_disease.csv",header=F)
#初始化两个矩阵
disease_third<-matrix(0,nrow = 169744,ncol = 1)
lncRNA_third<-matrix(0,nrow = 57600,ncol = 1)
Jaccard_disease<-matrix(as.matrix(Jaccard_disease))
Jaccard_lncRNA<-matrix(as.matrix(Jaccard_lncRNA))
GIP_disease<-matrix(as.matrix(GIP_disease))
GIP_lncRNA<-matrix(as.matrix(GIP_lncRNA))
#相同位置取最大值进行融合，并转换成数据框
for (i in 1:57600) {
  lncRNA_third[i]<-max(Jaccard_lncRNA[i,1],GIP_lncRNA[i,1])
}
lncRNA_third<-as.data.frame(matrix(lncRNA_third,nrow = 240,ncol = 240))
for (i in 1:169744) {
  disease_third[i]<-max(Jaccard_disease[i,1],GIP_disease[i,1])
}
disease_third<-as.data.frame(matrix(disease_third,nrow = 412,ncol = 412))








