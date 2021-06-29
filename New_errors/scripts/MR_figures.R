library(MendelianRandomization)

#Load file
#Plot the graph
#Save the graph

path<-"~/Desktop/Summer_project/Results/MR/server"
setwd(path)

results<-as.list(list.files(path))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

names=NULL
for (i in 1:length(results)){
  names[i]<-gsub(".RData.*","",results[[i]])
}

results_mr=NULL

for (i in 1:length(results)){
  if (i!=8){
    path<-"~/Desktop/Summer_project/Results/MR/server"
    setwd(path)
    mr<-loadRData(results[[i]])
    path<-"~/Desktop/Summer_project/Figures/MR"
    setwd(path)
    mr_plot(mr,interactive = FALSE)
    quartz.save(paste0(names[i],".png"),type="png")
    dev.off()
  } else {
    print("Error in lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,  : 
  NA/NaN/Inf in 'x'")
  }
}

mr.figures<-function(i){
  if (i!=8){
    path<-"~/Desktop/Summer_project/Results/MR/server"
    setwd(path)
    mr<-loadRData(results[[i]])
    path<-"~/Desktop/Summer_project/Figures/MR"
    setwd(path)
    mr_plot(mr,interactive = FALSE)
    quartz.save(paste0(names[i],".png"),type="png")
    dev.off()
  } else {
    print("Error in lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,  : 
  NA/NaN/Inf in 'x'")
  }
}

mr.figures(i=30)

path<-"~/Desktop/Summer_project/Results/MR/server"
setwd(path)
mr<-loadRData(results[[30]])
path<-"~/Desktop/Summer_project/Figures/MR"
setwd(path)
mr_plot(mr,interactive = FALSE)
quartz.save(paste0(names[30],".png"),type="png")
dev.off()


path<-"~/Desktop/Summer_project/Results/MR/server"
setwd(path)
mr<-loadRData(results[[2]])
path<-"~/Desktop/Summer_project/Figures/MR"
setwd(path)
mr_allmeth=mr_allmethods(mr)
mr_plot(mr_allmeth)
quartz.save(paste0(names[2],"_ALLMETH.png"),type="png")
dev.off()

x<-mr_ivw(mr)
x

mr<-loadRData(results[[6]])
