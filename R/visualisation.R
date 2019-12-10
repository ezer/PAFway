#' Draw network of enriched functional annotation pairs
#'
#' @param graph The output of either the pafway or pafway_edge_weight functions
#' @param pval The threshold of p-value at which to draw an arrow
#' @param adjMethod The method for correcting for multiple hypotheses.  This can be any method that is acceptable to the p.adjust function in the stats package: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none".  If this is NULL, then no adjustment will be made.
#' @param seed The random seed that will be used.
#' @return A matrix that has the same number of rows and columns as length(GOtypes).  This will contain p-values.
#' @export
#' @examples
#' a=matrix(c(0.1, 0.003, 0.005, 0.004, 0.5, 0.7, 0.001, 0.002, 0.003), nrow=3)
#' colnames(a)=c('A', 'B', 'C')
#' rownames(a)=c('A', 'B', 'C')
#' draw_network(a)
draw_network <- function(graph, pval = 0.05, adjMethod = NULL, seed = 123) {

    if (!is.null(adjMethod)) {
        graph2 = matrix(stats::p.adjust(graph, method = adjMethod), nrow = dim(graph)[1], ncol = dim(graph)[2])
        rownames(graph2) = rownames(graph)
        colnames(graph2) = colnames(graph)
        graph = graph2
    }

    GO.adj.matrix.norm.network <- graph
    GO.adj.matrix.norm.network[graph <= pval] <- 1
    GO.adj.matrix.norm.network[graph > pval] <- 0


    igraph.GO.allgenes <- igraph::graph.adjacency(GO.adj.matrix.norm.network, mode = "directed")
    igraph.GO.allgenes <- igraph::simplify(igraph.GO.allgenes)
    igraph::E(igraph.GO.allgenes)$col.path <- "grey"

    igraph::V(igraph.GO.allgenes)$size.degree <- igraph::degree(igraph.GO.allgenes)


    set.seed(seed)

    if(length(igraph::as_edgelist(igraph.GO.allgenes))==0){
        print("no edges are significant")
    }else{
    GGally::ggnet2(igraph::as_edgelist(igraph.GO.allgenes), label = TRUE, arrow.size = 7, arrow.gap = 0.025, node.size = "size.degree", edge.color = "grey50", layout.exp = 0.5) +
        ggplot2::ggtitle(paste0("Network of functional connectivity ")) + ggplot2::guides(size = FALSE)
    }
}
#' Draw network of enriched functional annotation pairs as a heatmap
#'
#' @param graph The output of either the pafway or pafway_edge_weight functions
#' @param adjMethod The method for correcting for multiple hypotheses.  This can be any method that is acceptable to the p.adjust function in the stats package: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none".  If this is NULL, then no adjustment will be made.
#' @param xlab The label for the x-axis of the heatmap
#' @param ylab The label for the y-axis of the heatmap
#' @param colPal The color palette of the heatmap
#' @return A matrix that has the same number of rows and columns as length(GOtypes).  This will contain p-values.
#' @export
#' @examples
#' nodes=paste("node", c(1:300))
#' set.seed(123)
#' randomGO=c("A", "B", "C", "D", "E", "F", "G", "H", "I",
#' "J", "K", "L", "M", "N")[sample(c(1:14), 300, replace=TRUE)]
#' names(randomGO)=nodes
#' edgesRandom=sapply(c(1:1000), function(i){
#'    nodes[sample(300, 2)]
#' })
#' getBinomPvalueRandom1=pafway(randomGO, t(edgesRandom), unique(randomGO))
#' draw_heatmap(getBinomPvalueRandom1)
#' colPal1=c(colorRampPalette(c("red3", "lightpink", "white", "white"))(20),
#' colorRampPalette(c("white", "white", "lightgreen", "darkgreen"))(20))
#' draw_heatmap(getBinomPvalueRandom1, adjMethod="bonferroni", xlab="Downstream",
#' ylab="Upstream", colPal=colPal1)
draw_heatmap <- function(graph, adjMethod = NULL, xlab = "downstream", ylab = "upstream", colPal = NULL) {
anythingSig=TRUE
    if (!is.null(adjMethod)) {
        graph2 = matrix(stats::p.adjust(graph, method = adjMethod), nrow = dim(graph)[1], ncol = dim(graph)[2])
        if(min(graph2)==1){
            print("nothing is even close to being significant")
            anythingSig=FALSE
        }
        rownames(graph2) = rownames(graph)
        colnames(graph2) = colnames(graph)
        graph = graph2
    }
if(anythingSig){
    if (is.null(colPal)) {
        gplots::heatmap.2(graph, trace = "none", key.title = "p-value", col = c(grDevices::colorRampPalette(c("red3", "lightpink",
            "white", "white"))(20), grDevices::colorRampPalette(c("white", "white", "lightblue", "blue4"))(20)), main = "",
            xlab = xlab, ylab = ylab, fig.height=10, fig.width=10, margins=c(15,15))
    } else {
        gplots::heatmap.2(graph, trace = "none", key.title = "p-value", col = colPal, main = "", xlab = xlab, ylab = ylab, margins=c(15,15))
    }
}
}

