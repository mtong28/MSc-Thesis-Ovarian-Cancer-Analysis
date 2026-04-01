gene <- read.csv(file = "./new.csv",header = TRUE, sep = ",")  # file= new.csv
gene$Time = as.numeric(gene$Time)

#install.packages("tidyverse")
#install.packages("rlang")

library("tidyverse")
library("survival")
## Define the function
#   l_LNC14112
FindBestCutoff_FromVectors <- function(Score, time, ind, PRINT=TRUE, 
                                       HighRisk = "High Risk",LowRisk="Low Risk"){
  Data <- cbind( Score = Score, Time = time, Ind = ind) %>% as.data.frame()
  Lower_q <- 0.2
  Upper_q <- 0.8
  
  Data_ordered <- Data[order(Data$Score),]
  
  q1 <- round(length(Data_ordered$Score)*Lower_q);
  q3 <- round(length(Data_ordered$Score)*Upper_q);
  
  surv <- Surv(Data_ordered$Time, Data_ordered$Ind)
  
  cutoffs <- seq(Lower_q,Upper_q, by=(Upper_q-Lower_q)/(q3-q1))
  
  pvals <- lapply(q1:q3, function(x)  {
    group <- ifelse(1:length(Data_ordered$Score) >= x , 1, 0)
    pval <- summary(coxph(surv ~ group))$sctest[3]
  }) %>% unlist()
  
  best_p <- min(pvals)
  best_cutoff <- cutoffs[pvals == best_p]
  best_int <- (q1:q3)[pvals == best_p]
  
  if (PRINT == TRUE) {
    print(cbind(cutoffs, -log2(pvals)) %>%
            data.frame() %>%
            setNames(c("cutoffs","p_values")) %>%
            ggplot(aes(x=cutoffs,y=p_values)) + 
            geom_point(pch="*", size=5) + 
            geom_vline(xintercept =  best_cutoff, linetype="dashed", colour="red") +
            geom_label(x=cutoffs[1], y=-log(best_p,2)/2, hjust=0,
                       label=paste("Best cut-off = ", round(best_cutoff,digits = 2),
                                   "\nLog-rank P = ", round(best_p,digits = 4),
                                   "\nInteger = ", best_int, sep = "")))
  }
  Data_ordered$BestGroup <- ifelse(1:length(Data_ordered$Score) >= best_int, HighRisk, LowRisk)
  Data_ordered[rownames(Data),]
}


### Visualization
group1 <- FindBestCutoff_FromVectors(gene$l_LNC14112,gene$Time,gene$Event)
group2 <- FindBestCutoff_FromVectors(gene$REV3L,gene$Time,gene$Event)
group3 <- FindBestCutoff_FromVectors(gene$l_LNC14112.REV3L,gene$Time,gene$Event)

#install.packages("survminer")
library("survminer")

# fit survival curves
require("survival")
# l_LNC14112
fit1 <- survfit(Surv(Time, Ind)~ BestGroup,data = group1 )
ggsurvplot(fit1, title = "l_LNC14112",font.title = c(35,"bold","darkblue"),
           legend = c(0.9,0.9),legend.title = "Group",legend.labs = c("High","Low"),
           pval = TRUE, break.time.by = 500,font.x=25, font.y=25,font.legend=30,pval.size=12,
           risk.table = TRUE, risk.table.col = "BestGroup", risk.table.y.text.col = TRUE)

# REV3L
fit2 <- survfit(Surv(Time, Ind)~ BestGroup,data = group2 )
ggsurvplot(fit2, title = "REV3L",font.title = c(35,"bold","darkblue"),
           legend = c(0.9,0.9),legend.title = "Group",legend.labs = c("High","Low"),
           pval = TRUE, break.time.by = 500,font.x=25, font.y=25,font.legend=30,pval.size=12,
           risk.table = TRUE, risk.table.col = "BestGroup", risk.table.y.text.col = TRUE)

# l_LNC14112/REV3L
fit3 <- survfit(Surv(Time, Ind)~ BestGroup,data = group3 )
ggsurvplot(fit3, title = "l_LNC14112/REV3L",font.title = c(35,"bold","darkblue"),
           legend = c(0.9,0.9),legend.title = "Group",legend.labs = c("High","Low"),
           pval = TRUE, break.time.by = 500,font.x=25, font.y=25,font.legend=30,pval.size=12,
           risk.table = TRUE, risk.table.col = "BestGroup", risk.table.y.text.col = TRUE)
