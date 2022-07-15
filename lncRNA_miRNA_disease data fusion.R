###############
#补全过程
#计算过程思想：例如lncRNA-miRNA与miRNA-disease乘积得到的矩阵，矩阵中的第[i,j]元素，如果lncRNA-miRNA的第i行与
#miRNA-disease的第j列的和都大于0，则第[i,j]除以lncRNA-miRNA的第i行与miRNA-disease的第j列的和,
#表示lncRNA通过miRNA与disease关联潜在可能性，这个过程类似于Jaccard相似系数的计算过程
#最后对得到的三个矩阵在第[i,j]取最大值进行融合，
#################
#加载数据集
lncRNA_disease<-read.csv(file = "D:/科研论文/实验材料/lncRNA_disease_240_412.csv",header=F)
lncRNA_name<-read.csv(file = "D:/科研论文/实验材料/lncRNAs_240.csv",header = F)
disease_name<-read.csv(file = "D:/科研论文/实验材料/diseases_412.csv",header = F)
lncRNA_miRNA<-read.csv(file = "D:/科研论文/实验材料/lncRNA_miRNA_240_495.csv",header=F)
miRNA_disease<-read.csv(file = "D:/科研论文/实验材料/miRNA_disease_495_412.csv",header=F)
#计算lncRNA-miRNA与miRNA-disease的乘积
lncRNA_disease_one<-as.data.frame(as.matrix(lncRNA_miRNA)%*%as.matrix(miRNA_disease))
#计算lncRNA-miRNA-disease的潜在关联矩阵
lncRNA_disease_one_pass<-matrix(0,nrow = 240,ncol = 412)#初始化一个240x412的0矩阵
row_sum<-as.matrix(rowSums(lncRNA_miRNA))#存储lncRNA-miRNA的行和
col_sum<-as.matrix(colSums(miRNA_disease))#存储lncRNA-miRNA的列和
##将计算得到的元素放到初始化矩阵的第[i，j]位置上
for (i in 1:240) {
  for (j in 1:412) {
  if((row_sum[i,1]>0)&(col_sum[j,1]>0)){
      lncRNA_disease_one_pass[i,j]<-lncRNA_disease_one[i,j]/(row_sum[i,1]+col_sum[j,1])
    }
  else{
      lncRNA_disease_one_pass[i,j]<-lncRNA_disease_one[i,j]
    }
  }
}
lncRNA_disease_pre<-matrix(as.matrix(lncRNA_disease))#先将原lncRNA-disease矩阵转换成一列
index<-which(lncRNA_disease_pre[,1]==0)#提取矩阵所有为0的元素的位置
lncRNA_disease_one_pre<-matrix(as.matrix(lncRNA_disease_one_pass))#先将原lncRNA_disease_one_pre矩阵转换成一列
#对lncRNA-miRNA矩阵的为0元素的位置用lncRNA_disease_one_pre相同位置的最大值进行几何填补
for (i in 1:96183) {
  lncRNA_disease_pre[index[i],1]<-lncRNA_disease_one_pre[index[i],1]
}
lncRNA_disease_pre<-as.data.frame(matrix(lncRNA_disease_pre,nrow = 240))









