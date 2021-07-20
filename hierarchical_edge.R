## Hierarchical edge bundling ##

rm(list=ls())

path<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

load("../Results/correlation_table_no_lewy.rds")
load("../Results/correlation_dataframe.rds")
load("../Results/pvals_no_lewy.rds")

df = data.frame(from = ifelse(gsub(".*_","",rownames(cor))!="lewy",gsub(".*_","",rownames(cor)),NA),
                to = ifelse(gsub("_.*","",rownames(cor))!="lewy",gsub("_.*","",rownames(cor)),NA),
                value = as.vector(cor$`cor[, 1]`),
                stringsAsFactors = FALSE)

df=df[complete.cases(df),]
df=cbind(df,p_val$p_val)

library(ggraph)

# plot
my_graph_hdl<-df[,c("from","to")]

vertices=data.frame(name = unique(c(as.character(df$from), as.character(df$to))),)

from <- match(my_graph_hdl$from, vertices$name)
to <- match(my_graph_hdl$to, vertices$name)

#labels
vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, df$from)))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves

my_graph<-graph_from_data_frame(my_graph_hdl,vertices=vertices)

ggraph(my_graph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point() +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", width=0.9,tension=1) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))



