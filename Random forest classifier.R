    #自动编码器特征选择
    disease_new<-cbind.data.frame(t(lncRNA_disease_pre),disease_third)
    lncRNA_new<-cbind.data.frame(lncRNA_disease_pre,lncRNA_third)
    network <-
      input() +
      dense(256, "sigmoid") +
      output("sigmoid")
    auto<-autoencoder(network, loss = "mean_squared_error")
    train_auto<-train(auto,as.matrix(lncRNA_new),validation_data=NULL,metrics=NULL,epochs=50,optimizer=keras::optimizer_rmsprop())
    learner_lncRNA<-encode(train_auto,as.matrix(lncRNA_new))
    network <-
      input() +
      dense(256, "sigmoid") +
      output("sigmoid")
    auto<-autoencoder(network, loss = "mean_squared_error")
    train_auto_disease<-train(auto,as.matrix(disease_new),validation_data=NULL,metrics=NULL,epochs=50,optimizer=keras::optimizer_rmsprop())
    learner_disease<-encode(train_auto_disease,as.matrix(disease_new))
    
      #数据的合并过程，disease，lncRNA的特征表示
      disease_new<-cbind.data.frame(disease_name,learner_disease)
      lncRNA_new<-cbind.data.frame(lncRNA_name,learner_lncRNA)
      lncRNA_disease_All<-merge(x=lncRNA_new,y=disease_new,by=NULL)
      
      #调节lncRNA和疾病列名的位置形成lncRNA-disease关联对
      d1<-subset(lncRNA_disease_All,select = 1)
      d2<-subset(lncRNA_disease_All,select = 256)
      d3<-subset(lncRNA_disease_All,select = c(-1,-258))
      lncRNA_disease_All<-data.frame(d1,d2,d3)
      
      #除去所有列和为0的列，，对实验来说没有任何影响
      d1 <- subset(lncRNA_disease_All,select=c(1,2))
      d2 <- subset(lncRNA_disease_All,select=c(-1,-2))
      d2 <- d2[,which(colSums(d2) > 0)] 
      lncRNA_disease_All <- cbind(d1, d2)
      
      #贴标签已经标签存在一个文件里面，每次用merge函数的结果lncRNA-disease pair是一样的
      LDExcl0_da <- read.csv("D:/科研论文/实验材料/全连接.csv",header = F)
      lncRNA_disease_new<-as.data.frame(cbind(LDExcl0_da[,1],lncRNA_disease_All))
      LDExcl0<-lncRNA_disease_new
      #更改前三列的列名，并进行归一化 
      colnames(LDExcl0)[1]<-"V1"
      colnames(LDExcl0)[2]<-"V2"
      colnames(LDExcl0)[3]<-"V3"
      B1=subset(LDExcl0[,], select=-c(V1,V2,V3))
      B2<-subset(LDExcl0[,], select=c(V1,V2,V3))
      normalize<-function(x){
        return((x-min(x))/(max(x)-min(x)))
      }
      B1<-as.data.frame(lapply(B1, normalize))
      LDExcl0 <- cbind(B2, B1)
      #构造正样本
      PositiveSample <- LDExcl0[LDExcl0[,1]==1, ]
      #构造无标签样本
      UnlabeledSample <- LDExcl0[LDExcl0[,1]==0, ]
      NB<-UnlabeledSample[-sp,c(-2,-3)]
      #构造负样本
      set.seed(12345)
      sp <- sample(nrow(UnlabeledSample), nrow(PositiveSample), replace = FALSE, prob = NULL)
      #sps<-sort(sp)
      NegativeSample <- UnlabeledSample[sp,]
      #构造正负样本矩阵 
      TrainingSample <- rbind(PositiveSample, NegativeSample)
      #xgboost算法——十折交叉验证
      PTB<-TrainingSample[1:2695,c(-2,-3)]
      NTB<-TrainingSample[2698:5392,c(-2,-3)]    
      #从样本集中每次随机抽取269个正负样本合成训练集，剩下的作为训练集
      #sample_number生成2697个0的列表
      #sample_num初试化测试集的长度，
      sample_number<-matrix(0,nrow = 2695,ncol = 1)
      sample_num<-matrix(0,nrow = 539,ncol = 1)
      for (i in 1:5) {
        sample_set<-sample(which(sample_number[,1]==0),539,replace = F)#生成随机数
        sample_number[sample_set,1]<-i#在sample_number中随机数的位置换成循环的次数i
        sample_num<-cbind.data.frame(sample_num,sample_set)#将随机数合并到一起，第九次之后剩下的所有0 的位置作为第十次的测试集的个数
      }
      sumauc <- 0#auc初试化
      i<-1
      sample_PT<-PTB[sample_num[,i+1],]#将sample_number剩余0的位置，在PTB中取同样的位置作为正样本测试集
      sample_NT<-NTB[sample_num[,i+1],]#将sample_number剩余0的位置，在NTB中取同样的位置作为负样本测试集
      sample_test<-rbind.data.frame(sample_PT,sample_NT,NB)#合并成测试集
      value<-sample_test$V1#测试集的标签
      for(i in 1:5)
      {
        sample_PT<-PTB[sample_num[,i+1],]#将sample_number剩余0的位置，在PTB中取同样的位置作为正样本测试集
        sample_NT<-NTB[sample_num[,i+1],]#将sample_number剩余0的位置，在NTB中取同样的位置作为负样本测试集
        sample_test<-rbind.data.frame(sample_PT,sample_NT,NB)#合并成测试集
        # value<-sample_test$V1#测试集的标签
        sample_train<-rbind.data.frame(PTB[-sample_num[,i+1],],NTB[-sample_num[,i+1],])#构造训练集
        #算法randomforest
        rf=randomForest(as.factor(V1)~.,data = sample_train, mtry=floor(sqrt(512)), importance = TRUE, ntree=500, proximity=TRUE)
        pred <- predict(rf, sample_test[,-1],type="prob")
        pred<-as.data.frame(pred)
        value<-cbind.data.frame(value,pred[,2])
        pred1 <- prediction(pred[,2], sample_test$V1)
        roc <- performance(pred1, "tpr", "fpr")
        auc <- performance(pred1, "auc")@y.values
        sumauc <- sumauc + as.numeric(auc[[1]])
      }
      print(sumauc/5)



