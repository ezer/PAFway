#### Public ones
#' Find pairwise-associations between annotations in a network without edge weights.
#'
#' @param GO A vector of Strings, equal to the length of all the nodes.  The names of the vector should be the names of the nodes.  The values should either be the functional annotations, or a concatenation of the functional annotations, separated by a "_" symbol.
#' @param edges A matrix of Strings, with at least two columns.  Each row will represent an edge, linking the node in the first column to the node in the second column.  Please make sure the node names are the same as those in "GO"
#' @param GOtypes This is a vector that contains the functional annotations or GO terms that are of interest
#' @param exact A boolean.  If it is true, it will look for an exact match between the term in GOtypes and the vector GO.  Otherwise, it will look for substrings.
#' @param adjustByEdgeCount A boolean.  If true, then the probability of observing a functional annotation will be calculated in terms of the number of nodes, but if it is false, it is calculated in terms of the number of edges that contain that node.
#' @return A matrix that has the same number of rows and columns as length(GOtypes).  This will contain p-values.
#' @export
#' @examples
#' nodes=paste("node", c(1:10))
#' set.seed(123)
#' randomGO=c("A", "B", "C")[sample(c(1:3), 10, replace=TRUE)]
#' names(randomGO)=nodes
#' edgesRandom=sapply(c(1:100), function(i){
#'    nodes[sample(10, 2)]
#' })
#' pafway(randomGO, t(edgesRandom), unique(randomGO))
pafway <- function(GO, edges, GOtypes, exact = TRUE, adjustByEdgeCount = FALSE) {

    GOinNetwork = GO[unique(c(edges[, 1], edges[, 2]))]
    sapply(GOtypes, function(i) {
        sapply(GOtypes, function(j) {

            if (length(which(grepl(i, GOinNetwork[edges[, 1]]) == 0)) & length(which(grepl(i, GOinNetwork[edges[,
                2]]) == 0))) {

            }


            if (exact) {
                a = length(which(GO[edges[, 1]] == i & GO[edges[, 2]] == j))


                if (adjustByEdgeCount) {
                  p_bot = (length(which(GOinNetwork[edges[, 1]] == i))/length(edges[, 1])) * (length(which(GOinNetwork[edges[,
                    2]] == j))/length(edges[, 2]))

                } else {
                  p_bot = (length(which(GOinNetwork == i))/length(GOinNetwork)) * (length(which(GOinNetwork ==
                    j))/length(GOinNetwork))
                   }

            } else {
                a = length(which(grepl(i, GOinNetwork[edges[, 1]]) & grepl(j, GOinNetwork[edges[, 2]])))

                if (adjustByEdgeCount) {
                  p_bot = (length(grep(i, GOinNetwork[edges[, 1]]))/length(edges[, 1])) * (length(grep(j, GOinNetwork[edges[,
                    2]]))/length(edges[, 2]))
                } else {
                  p_bot = (length(grep(i, GOinNetwork))/length(GOinNetwork)) * (length(grep(j, GOinNetwork))/length(GOinNetwork))
                }
              }

            b = stats::binom.test(a, length(edges[, 1]), p = p_bot, alternative = c("greater"), conf.level = 0.95)
            b$p.value
        })
    })
}

#' Find pairwise-associations between annotations in a network with edge weights.
#'
#' @param GO A vector of Strings, equal to the length of all the nodes.  The names of the vector should be the names of the nodes.  The values should either be the functional annotations, or a concatenation of the functional annotations, separated by a "_" symbol.
#' @param edges A matrix of Strings, with at least three columns.  Each row will represent an edge, linking the node in the first column to the node in the second column, and the third column will contain an edge weight.  Please make sure the node names are the same as those in "GO"
#' @param GOtypes This is a vector that contains the functional annotations or GO terms that are of interest
#' @param exact A boolean.  If it is true, it will look for an exact match between the term in GOtypes and the vector GO.  Otherwise, it will look for substrings.
#' @param adjustByEdgeCount A boolean.  If true, then the probability of observing a functional annotation will be calculated in terms of the number of nodes, but if it is false, it is calculated in terms of the number of edges that contain that node.
#' @param step FFT is used to speed up the calculation.  In the first step, a certain number of values are evenly sampled from the function, across its range.  This value will determine the distance between sampled points.
#' @param thresholdZero In order to decrease the space and time requirements, values in the probability distributions that are below a certain threshold are set to be exactly zero.  This is the threshold.
#' @return A matrix that has the same number of rows and columns as length(GOtypes).  This will contain p-values.
#' @examples
#' nodes=paste("node", c(1:10))
#' set.seed(123)
#' randomGO=c("A", "B", "C")[sample(c(1:3), 10, replace=TRUE)]
#' names(randomGO)=nodes
#' edgesRandom=sapply(c(1:20), function(i){
#'    nodes[sample(10, 2)]
#' })
#' getBinomPvalueRandom1=pafway_edge_weights(randomGO, cbind(t(edgesRandom),
#' rnorm(length(edgesRandom[1,]), 1, 0.001)), unique(randomGO))
#' @export
pafway_edge_weights <- function(GO, edges, GOtypes, exact = TRUE, adjustByEdgeCount = FALSE, step = 0.001, thresholdZero = 1e-04) {

    # gosums=
    params = lapply(GOtypes, function(i) {
        sapply(GOtypes, function(j) {
            if (exact) {
                z = length(which(GO[edges[, 1]] == i & GO[edges[, 2]] == j))
                if (adjustByEdgeCount) {
                  p_bot = (length(which(GO[edges[, 1]] == i))/length(edges[, 1])) * (length(which(GO[edges[,
                    2]] == j))/length(edges[, 2]))

                } else {
                  p_bot = (length(which(GO == i))/length(GO)) * (length(which(GO == j))/length(GO))
                  }
            } else {
                z = length(which(grepl(i, GO[edges[, 1]]) & grepl(j, GO[edges[, 2]])))
                if (adjustByEdgeCount) {
                  p_bot = (length(grep(i, GO[edges[, 1]]))/length(edges[, 1])) * (length(grep(j, GO[edges[, 2]]))/length(edges[,
                    2]))
                } else {
                  p_bot = (length(grep(i, GO))/length(GO)) * (length(grep(j, GO))/length(GO))
                }
            }

            c(p_bot, z)
        })
    })
    param_vector = do.call(cbind, params)
    den = stats::density(as.numeric(edges[, 3]))
    minDen = min(den$x)

    maxDen = max(den$x)

    denFun = stats::approxfun(stats::density(as.numeric(edges[, 3])), yleft = 0, yright = 0)

    cdfFun = stats::ecdf(as.numeric(edges[, 3]))

    range = seq(minDen - step * 10, maxDen + 10 * step, step)  #make this a parameter
    denSamp = sapply(range, function(i) {
        denFun(i)
    })
    cdfSamp = sapply(range, function(i) {
        cdfFun(i)
    })
    denSampOr = denSamp
    rangeOr = range

    slides = min(range)
    sum = rep(0, length(param_vector[1, ]))

    for (k in c(1:length(edges[, 3]))) {

        probKedges = sapply(param_vector[1, ], function(itemp) {
            stats::dbinom(k, length(edges[, 3]), itemp)
        })

        pZgXk = stats::convolve(denSamp, rev(cdfSamp), type = "open")
        pZgXk = pZgXk/max(pZgXk)
        rangeExpanded = c(range, seq(max(range) + step, max(range) + step * (length(range) - 1) + 1e-10, step)) +
            slides

        ## chop off artifact
        firstOne = which(pZgXk == 1)[1]
        sum = sum + sapply(c(1:length(param_vector[2, ])), function(itemp_id) {
            probKedges[itemp_id] * stats::approxfun(rangeExpanded[1:firstOne], pZgXk[1:firstOne], yleft = 0, yright = 1)(param_vector[2,
                itemp_id])
        })

        denSamp = stats::convolve(denSampOr, rev(denSamp), type = "open")
        denSamp = denSamp/max(denSamp)

        # update range, denSamp, cdfSamp
        minDen = rangeExpanded[min(which(denSamp > thresholdZero))]
        maxDen = rangeExpanded[max(which(denSamp > thresholdZero))]
        range = seq(minDen - step * 10, maxDen + 10 * step, step)  #make this a parameter
        denFun2 = stats::approxfun(rangeExpanded, denSamp, yleft = 0, yright = 0)
        denSamp = sapply(range, function(i) {
            denFun2(i)
        })

        if ((length(range) - length(denSampOr)) < 0) {
            denSampOr = denSampOr[1:length(range)]
        } else {
            denSampOr = c(denSampOr, rep(0, length(range) - length(denSampOr)))
        }

        if ((length(range) - length(cdfSamp)) < 0) {
            cdfSamp = cdfSamp[1:length(range)]
        } else {
            cdfSamp = c(cdfSamp, rep(1, length(range) - length(cdfSamp)))
        }

    }

    sum = 1 - sum


    gosums = matrix(sum, nrow = sqrt(length(sum)), ncol = sqrt(length(sum)))
    rownames(gosums) = GOtypes
    colnames(gosums) = GOtypes
    return(gosums)
}




# #### Private ones
#
# make_pafway_object <- function(edgelist.filename, GO.filename, genelist.filename = NULL, threshold.edgeImportance = 0.022,
#     sep = ",", header = TRUE) {
#     temp1 = parse_network(edgelist.filename, threshold.edgeImportance, sep = sep, header = header)
#
#     if (is.null(genelist.filename)) {
#         genelist = unique(c(unique(temp1[, 1]), unique(temp1[, 2])))
#     } else {
#         genelist <- read.csv(genelist.filename, header = TRUE, stringsAsFactors = FALSE, sep = sep)
#     }
#     temp2 = parse_GO_file(GO.filename, genelist = genelist)
#     if (length(which(!(temp2[, 1] %in% c(temp1[, 1], temp1[, 2])))) != 0) {
#         temp2 = temp2emp[which(temp2[, 1] %in% c(temp1[, 1], temp1[, 2])), ]
#     }
#     return(list(network = temp1, description = temp2))
# }
#
# parse_network <- function(edgelist.filename, threshold.edgeImportance = 0.022, sep = ",", header = TRUE) {
#     res.edgeImportance <- read.csv(edgelist.filename, header = header, stringsAsFactors = FALSE, sep = sep)
#
#     ## if no edge weights, include all####@TODO
#
#     ## -- exclude edges below cutoff
#     res.edgeImportance.filt <- res.edgeImportance[res.edgeImportance[, 3] > threshold.edgeImportance, ]
#     ## -- IGRAPH object vertices.available.edgeImportance <- unique(c(res.edgeImportance.filt[,1],
#     ## res.edgeImportance.filt[,2])) igraph.edgeImportance <- graph.data.frame(res.edgeImportance.filt,
#     ## vertices.available.edgeImportance, directed=T) igraph.edgeImportance <-
#     ## induced_subgraph(igraph.edgeImportance, components(igraph.edgeImportance)$membership==1) edge.list.GO.link
#     ## <- igraph::get.edgelist(igraph.edgeImportance) print('dimensions') print(dim(edge.list.GO.link))
#     ## print(length(res.edgeImportance.filt[,3]))
#     return(res.edgeImportance.filt)
# }
#
# parse_GO_file <- function(GO.filename, genelist, sep = "\t", messy = T, colGeneName = 1, colAnnotation = 5) {
#     if (!messy) {
#         gene.list.GO <- read.csv(GO.filename, skip = 1, stringsAsFactors = FALSE, sep = sep, header = FALSE)
#     } else {
#         gene.list.GO <- read.csv(GO.filename, skip = 1, stringsAsFactors = FALSE, sep = sep, header = FALSE)[,
#             c(1, 5)]
#         gene.list.GO <- t(sapply(genelist, function(i) {
#             c(i, paste(gene.list.GO[which(gene.list.GO[, 1] == i), 2], collapse = "_"))
#         }))
#         colnames(gene.list.GO) = c("gene", "GO")
#         return(gene.list.GO)
#     }
#     important.gene.list.GO = gene.list.GO[which(gene.list.GO[, 1] %in% gene.list[[1]]), ]
#
#     gene.list.GO.modified <- matrix(NA, length(important.gene.list.GO[, 1]), 2, dimnames = list(NULL, c("gene",
#         "GO")))
#     gene.list.GO.modified[, "gene"] <- important.gene.list.GO[, 1]
#     gene.list.GO.modified[, "GO"] <- apply(important.gene.list.GO[, 2:(dim(important.gene.list.GO)[2])], 1, function(x) {
#         x <- paste(x, collapse = "_")
#         x <- unique(unlist(strsplit(x, split = "_")))
#         paste(x, collapse = "_")
#     })
#
#     return(gene.list.GO.modified)
# }
#
# ####Potentially delete:
#
# make_pafway_object_file <- function(edgelist.filename, GO.filename, genelist.filename = NULL, GO.terms.interest = NA,
#                                     numPert = 200, threshold.edgeImportance = 0.022, sep = ",", header = TRUE) {
#     obj = make_pafway_object(edgelist.filename = edgelist.filename, GO.filename = GO.filename, genelist.filename = genelist.filename,
#                              threshold.edgeImportance = threshold.edgeImportance, sep = sep, header = header)
#     print(dim(obj[["network"]]))
#     return(obj)
#     # return(enrichNet(data=obj, GO.terms.interest=GO.terms.interest, numPert=numPert))
# }
#
# filter_nodes_wo_functions <- function(data, GO.terms.interest, exact = F) {
#     # @TODO exact
#     genesWithGO = data[["description"]][which(apply(sapply(GO.terms.interest, function(j) {
#         grepl(j, data[["description"]][, 2])
#     }), 1, function(i) {
#         sum(i) > 0
#     })), 1]
#     subEdges = data[["network"]][which((data[["network"]][, 1] %in% genesWithGO) & (data[["network"]][, 2] %in%
#                                                                                         genesWithGO)), ]
# }
#
