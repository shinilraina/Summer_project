## The files you need for this script to run: 
# 1. The HDL correlation table
# 2. The matrix of 'n' values for the pc algorithm
# 3. The file of MR results

# Packages: pcalg, ggraph, igraph

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

pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(correlation1),alpha = 0.05)
#pc.fit1 <- pc(suffStat, indepTest = gaussCItest, p = ncol(correlation1),alpha = 0.05/(ncol(correlation1)*nrow(correlation1)))

labels<-c("ALS","AD","E","MS","PD","S")
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

#set.seed(123)
#a<-as.directed(graph_from_graphnel(pc.fit@graph))
#labels<-c("ALS","AD","E","MS","PD","S")
#V(a)$size<-35
#V(a)$labels<-labels
#V(a)$color<-"white"
#E(a)$arrow.size<-0.5
#E(a)$arrow.width<-0.75

library(ggraph)
edges<-as.matrix(showEdgeList(pc.fit))
edges<-edges[2]

edge_list<-matrix(nrow=12,ncol=2)
edge_list[,1]<-edges[[1]][1:12] #from
edge_list[,2]<-edges[[1]][13:length(edges[[1]])] #to
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

ida_list<-merge(as.data.frame(edge_list),ida_df,by=c("from","to"))
effect_size<-ida_list$ida_effect

#Curved/final plot for HDL-CI plot
labels<-c("ALS","S","PD","E","AD","MS")
ci_plot<-ggraph(my_graph, layout = 'linear', circular = TRUE) +
  geom_edge_hive(aes(edge_colour=effect_size,start_cap = circle(5, 'mm'),
                     end_cap = circle(5, 'mm')),alpha=1,
                 arrow = arrow(length = unit(0.8, 'mm'),type = "closed"),lineend = "round",strength=0.5)+
  geom_node_point(aes(x=x*1, y=y*1),size=10,colour="grey",alpha=0.005)+
  geom_node_text(aes(label=labels,x=x*1.0,y=y*1.0), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2) +
  theme_void()+
  scale_edge_colour_distiller(palette ="RdBu")+
  labs(edge_color="Effect size")+
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggsave(ci_plot,filename="../Figures/CI/ci_hdl_network_ggraph.png")

## Undirected HDL-CI plot using edge_arc
labels<-c("ALS","AD","E","MS","S","PD")
undirected<-ggraph(my_graph, layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(edge_width=10*abs(effect_size),alpha=0.4,edge_color=effect_size))+
  guides(edge_alpha = "none", edge_width = "none")+
  geom_node_point(aes(x=x*0.95, y=y*0.95),size=2,colour="white",alpha=1)+
  geom_node_text(aes(label=labels,x=x*1.09,y=y*1.09), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2) +
  theme_void()+
  scale_edge_colour_distiller(palette ="RdBu")+
  labs(edge_color="Effect size")+
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave(undirected,filename="../Figures/CI/ci_hdl_undirected_ggraph.png")

#Curved/final MR plot
load("../Results/MR/MR_results.rds")

mr_graph<-mr_table[mr_table$"Pvalue"<0.05,c("Exposure","Outcome","Estimate","Pvalue")]
colnames(mr_graph)<-c("from","to","effect","pval")

mr_plot<-graph_from_data_frame(mr_graph)
plot(mr_plot,layout=layout.circle(mr_plot))

label_mr<-c("ALS","AD","E","LBD","PD","S","MS")

ggraph(mr_plot, layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(edge_width=mr_graph$effect,edge_colour=mr_graph$effect,start_cap = circle(5, 'mm'),
                     end_cap = circle(5, 'mm')),alpha=0.5,
                 arrow = arrow(length = unit(mr_graph$effect, 'mm'),type = "closed"),lineend = "round")+
  geom_node_point(aes(x=x*1, y=y*1),size=10,colour="grey",alpha=0.005)+
  geom_node_text(aes(label=label_mr,x=x*1.0,y=y*1.0), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2) +
  theme_void()+
  scale_edge_colour_distiller(palette ="RdBu")+
  labs(edge_color="Effect size")+
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggsave(ci_plot,filename="../Figures/CI/MR_network_ggraph.png")

#Combine MR and HDL graphs together
mr_graph$group<-rep("MR",times=nrow(mr_graph))
mr_list<-mr_graph[,-4]
ida_list$group<-rep("HDL",times=nrow(ida_list))
colnames(ida_list)<-c("from","to","effect","group")

mr_list$from<-ifelse(mr_list$from=="alz","AD",
                ifelse(mr_list$from=="als","ALS",
                       ifelse(mr_list$from=="lewy","LBD",
                              ifelse(mr_list$from=="stroke","S",
                                     ifelse(mr_list$from=="p","PD",
                                            ifelse(mr_list$from=="ms","MS",
                                                   ifelse(mr_list$from=="epilepsy","E",NA)))))))
mr_list$to<-ifelse(mr_list$to=="alz","AD",
              ifelse(mr_list$to=="als","ALS",
                     ifelse(mr_list$to=="lewy","LBD",
                            ifelse(mr_list$to=="stroke","S",
                                   ifelse(mr_list$to=="p","PD",
                                          ifelse(mr_list$to=="ms","MS",
                                                 ifelse(mr_list$to=="epilepsy","E",NA)))))))

comparison<-rbind(ida_list,mr_list)

mr_plot<-graph_from_data_frame(comparison)
true<-ifelse(comparison$group=="HDL","bold","bold")

label_mr<-c("AD","ALS","MS","PD","S","LBD","E")

comparison_plot<-ggraph(mr_plot, layout = 'linear', circular = TRUE) +
  geom_edge_link(aes(start_cap = circle(5, 'mm'),
                     end_cap = circle(5, 'mm'),
                     colour=ifelse(comparison$group=="MR","blue","red")),alpha=0.9,
                 arrow = arrow(length = unit(1.8, 'mm'),type = "closed"),
                 lineend = "round",show.legend=FALSE)+
  geom_node_point(aes(x=x*1, y=y*1),size=10,colour="grey",alpha=0.005)+
  geom_node_text(aes(label=label_mr,x=x*1.0,y=y*1.0), 
                 size=3.5, alpha=0.8,repel=FALSE,fontface = 2) +
  theme(legend.direction = "vertical",legend.position = "left",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  theme_void()

ggsave(comparison_plot,filename="../Figures/CI/comparison_ggraph.png")
