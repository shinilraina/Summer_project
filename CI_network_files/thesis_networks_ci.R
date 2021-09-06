# The files you need for this script to run: 
# 1. The HDL correlation table
# 2. The matrix of 'n' values for the pc algorithm
# 3. The file of MR results

# Packages: pcalg, ggraph, igraph, cowplot

rm(list=ls())

library(pcalg)

path<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

load("~/Desktop/Summer_project/Results/correlation_table_no_lewy.rds")

num_matrix<-read.csv("../Data/num_trial.csv",header=TRUE)
rownames(num_matrix)<-num_matrix[,1]
num_matrix<-num_matrix[,-1]

pc_cors<-as.matrix(correlation1)

suffStat<-list(C=pc_cors,n=min(as.matrix(num_matrix),na.rm=TRUE))

pc.fit1 <- pc(suffStat, indepTest = gaussCItest, p = ncol(correlation1),alpha = 0.05,verbose=TRUE)
n=sum(pc.fit1@n.edgetests)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(correlation1),alpha = 0.05/n,verbose=TRUE)

labels<-c("ALS","AD","E","MS","PD","S")
a<-function(a){jointIda(a,vars, as.matrix(correlation1), pc.fit@graph)}
b<-function(a){idaFast(a,vars, as.matrix(correlation1), pc.fit@graph)}
vars<-c(1,2,3,4,5,6)

ida_results<-lapply(vars, b)
ida<-data.frame()
for (i in 1:length(ida_results)){
  for (j in 1:length(ida_results)){
    ida[i,j]=ida_results[[i]][j]
  }
}
colnames(ida)<-labels
rownames(ida)<-labels

ida_df<-stack(ida[1:ncol(ida)])
colnames(ida_df)<-c("ida_effect","to")
ida_df$from<-rep(labels,times=6)
ida_df<-ida_df[which(ida_df$ida_effect!=1 & ida_df$ida_effect!=0),]

library(graph)
library(igraph)

library(ggraph)
edges<-as.matrix(showEdgeList(pc.fit))
edges[1] #Undirected edge 4-5
edges<-edges[2]

edge_list<-matrix(nrow=13,ncol=2)
edge_list[12,c(1,2)]<-c(4,5)
edge_list[13,c(1,2)]<-c(5,4)
edge_list[c(1:11),1]<-edges[[1]][1:11] #from
edge_list[c(1:11),2]<-edges[[1]][12:length(edges[[1]])] #to

colnames(edge_list)<-c("from","to")

edge_list[,"from"]<-ifelse(edge_list[,"from"]==1,"ALS",
                           ifelse(edge_list[,"from"]==2,"AD",
                                  ifelse(edge_list[,"from"]==3,"E",
                                         ifelse(edge_list[,"from"]==4,"MS",
                                                ifelse(edge_list[,"from"]==5,"PD",
                                                       ifelse(edge_list[,"from"]==6,"S",NA))))))

edge_list[,"to"]<-ifelse(edge_list[,"to"]==1,"ALS",
                         ifelse(edge_list[,"to"]==2,"AD",
                                ifelse(edge_list[,"to"]==3,"E",
                                       ifelse(edge_list[,"to"]==4,"MS",
                                              ifelse(edge_list[,"to"]==5,"PD",
                                                     ifelse(edge_list[,"to"]==6,"S",NA))))))

my_graph<-graph_from_edgelist(edge_list)

ida_df<-ida_df[,c(3,2,1)]
ida_list<-merge(ida_df,as.data.frame(edge_list),by=c("from","to"))
effect_size<-ida_list$ida_effect

## Directed HDL-CI plot using edge_arc
set.seed(123)
directed<-ggraph(my_graph, layout = 'linear', circular = TRUE) +
  geom_node_point(aes( x = x*1.05, y=y*1.05)) +
  scale_edge_colour_distiller(palette ="RdBu")+
  geom_edge_arc(aes(edge_width=abs((effect_size)), 
                    edge_color=(effect_size)),
                arrow = arrow(length = unit(c(0.1,0.1,0.1) ,"inches"),type = "closed"),
                alpha=0.3,lineend = "round",end_cap=circle(2,"mm"),start_cap=circle(2,"mm"))+
  geom_edge_arc(aes(edge_color=(effect_size),edge_width=0.0001),
                arrow = arrow(length=unit(0.7,"mm"),
                              type = "closed"),
                alpha=1,end_cap=circle(20,"mm"),start_cap=circle(20,"mm"))+
  guides(edge_alpha = "none", edge_width = "none")+
  geom_node_point(aes(x=x*0.99, y=y*0.99),size=2,colour="black",alpha=0)+
  geom_node_text(aes(label=name,x=x*1.15,y=y*1.15), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2) +
  theme_void()+
  labs(edge_color="Effect size")+
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggsave(directed,filename="../Figures/CI/ci_hdl_directed_ggraph.png")

#Curved/final MR plot
load("../Results/MR/MR_results.rds")

mr_graph<-mr_table[mr_table$"Pvalue"<0.05,c("Exposure","Outcome","Estimate","Pvalue","adj.sig")]
colnames(mr_graph)<-c("from","to","effect","pval","adj")

mr_graph[,"from"]<-ifelse(mr_graph[,"from"]=="als","ALS",
                          ifelse(mr_graph[,"from"]=="alz","AD",
                                 ifelse(mr_graph[,"from"]=="epilepsy","E",
                                        ifelse(mr_graph[,"from"]=="ms","MS",
                                               ifelse(mr_graph[,"from"]=="p","PD",
                                                      ifelse(mr_graph[,"from"]=="stroke","S",
                                                             ifelse(mr_graph[,"from"]=="lewy","LBD",NA)))))))

mr_graph[,"to"]<-ifelse(mr_graph[,"to"]=="als","ALS",
                        ifelse(mr_graph[,"to"]=="alz","AD",
                               ifelse(mr_graph[,"to"]=="epilepsy","E",
                                      ifelse(mr_graph[,"to"]=="ms","MS",
                                             ifelse(mr_graph[,"to"]=="p","PD",
                                                    ifelse(mr_graph[,"to"]=="stroke","S",
                                                           ifelse(mr_graph[,"to"]=="lewy","LBD",NA)))))))

mr_graph<-mr_graph[c(2,1,3,7,4,5,6),]
mr_plot<-graph_from_data_frame(mr_graph)

#MR only arc graph
width<-ifelse(mr_graph$adj=="Significant",abs(mr_graph$effect),0.0001)
colours<-ifelse(mr_graph$adj=="Significant",exp(mr_graph$effect),0.001)
alphas<-ifelse(mr_graph$adj=="Significant",0.4,0.001)
mr_graph$linetype=ifelse(mr_graph$adj=="Significant","bold","dashed")

mr<-ggraph(mr_plot, layout = 'linear', circular = TRUE) +
  geom_node_point(aes( x = x*1.05, y=y*1.05))+
  geom_edge_arc(aes(edge_width=width,
                    colour=colours,linetype=mr_graph$linetype,alpha=alphas),
                arrow = arrow(length = unit(c(0.1,0.1,0.1) ,"inches"), 
                              type = "closed"),
                lineend = "round",end_cap=circle(2,"mm"),start_cap=circle(2,"mm"))+
  guides(edge_alpha = "none", edge_width = "none")+
  geom_node_point(aes(x=x*0.99, y=y*0.99),size=2,colour="black",alpha=0)+
  geom_node_text(aes(label=name,x=x*1.15,y=y*1.15), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2)+
  scale_edge_colour_distiller(palette ="RdBu")+
  theme_void() +
  labs(edge_color="Odds ratio")+
  scale_edge_linetype(name="Significance",
                      labels=c("BH adjusted significant","Nominally significant"))+
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

library(cowplot)
pdf("../Figures/CI/comparison_side_ggraph.pdf", width = 6, height = 9)
plot_grid(directed, mr, labels = "AUTO",nrow=2,ncol=1)
dev.off()