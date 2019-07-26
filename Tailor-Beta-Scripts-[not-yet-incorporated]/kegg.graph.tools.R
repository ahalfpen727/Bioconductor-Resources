#-------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#-------------------------------------------------------------------------------
displayNetwork <- function(tbl.kegg, layout.file, gene.type.map,
                           genes.of.interest=NULL)
{
    if (!is.null(genes.of.interest)) {
        from.rows <- which(tbl.kegg$from %in% genes.of.interest)
        to.rows   <- which(tbl.kegg$to %in% genes.of.interest)
        keepers <- intersect(from.rows, to.rows)
        stopifnot(length(keepers) > 0)
        tbl.kegg <- tbl.kegg[keepers,]
        }

    kegg.geneIDs <- unique(c(tbl.kegg$from, tbl.kegg$to))
    geneID.geneSymbol.map <-
       AnnotationDbi::mget(kegg.geneIDs, envir=org.Hs.egSYMBOL, ifnotfound=NA)
    g <- createKeggGraph(tbl.kegg, geneID.geneSymbol.map, gene.type.map)
    cw <- new.CytoscapeWindow("KEGG Cancer-Related", g, overwriteWindow=TRUE)
    displayGraph(cw)
    hideAllPanels(cw)
    restoreLayout(cw, layout.file)
    fitContent(cw)
    setZoom(cw, 0.8 * getZoom(cw))
    cw

} # displayNetwork
#-------------------------------------------------------------------------------
combineKeggPathways <- function(pathways, organism)
{
    tbl.net <- data.frame(stringsAsFactors=FALSE)
    for (pathway in as.character(pathways)) {
        filename <- sprintf("%s-%s.kgml", organism, pathway)
        if (!filename %in% list.files())
           retrieveKGML(pathway, organism=organism, destfile=filename,
                         method="internal")
        tbl.new <- parseKGML2DataFrame(filename)
        tbl.new$from <- as.character(tbl.new$from)
        tbl.new$to <- as.character(tbl.new$to)
        tbl.new$subtype <- as.character(tbl.new$subtype)
        tbl.new$from <- sapply(strsplit(tbl.new$from, ":"),
                              function(tokens) return(tokens [2]))
        tbl.new$to <- sapply(strsplit(tbl.new$to, ":"),
                            function(tokens) return(tokens [2]))
        pathway.name <- sprintf("kegg:%s", pathway)
        tbl.new$source <- rep(pathway.name, nrow(tbl.new))
        tbl.net <- unique(rbind(tbl.net, tbl.new))
        } # for pathway

  tbl.net

} # combineKeggPathways
#-------------------------------------------------------------------------------
tumorViz <- function(cw, tumor.name=NULL)
{
    setNodeTooltipRule(cw, "nodeType")
    setNodeSizeRule(cw, "score", c(0, 1, 3, 5, 20), c(50, 80, 120, 160, 300),
                    mode="interpolate")
    setNodeColorRule(cw, "lfc",  c(-3,0,3),  c("#00FF00", "#FFFFFF", "#FF0000"),
                     mode="interpolate", default.color="#0000FF")
    setNodeLabelRule(cw, "label");
    types <- c("gene","growthFactor","info","kinase","kinase TF","TF",
              "receptor kinase")
    shapes <- c("ellipse", "vee", "ellipse", "hexagon", "rect", "rect", "diamond")
    setNodeShapeRule(cw, "nodeType", types, shapes); redraw(cw)

    setNodeBorderColorRule(cw, "copyNumber", c("-2", "-1", "0", "1", "2"),
                           c("#000000", "#000000", "#000000", "#0000FF", "#0000FF"),
                           mode="lookup")
    setNodeBorderWidthRule(cw, "copyNumber", c("-2", "-1", "0", "1", "2"),
                           c(40, 20, 1, 20, 40))

    g <- cw@graph

    edge.map <-  list(
      "activation"                 = list(line="SOLID",          targetArrow="Arrow",    sourceArrow="No Arrow", color="#00AA00", width=6),
      "activation;phosphorylation" = list(line="SINEWAVE",       targetArrow="Arrow",    sourceArrow="No Arrow", color="#00AA00", width=3),
      "binding/association"        = list(line="SOLID",          targetArrow="Arrow",    sourceArrow="Arrow",    color="#00AAAA", width=3),
      "inhibition;phosphorylation" = list(line="SINEWAVE",       targetArrow="T",        sourceArrow="No Arrow", color="#DD0000", width=3),
      "inhibition"                 = list(line="DASH_DOT",       targetArrow="T",        sourceArrow="No Arrow", color="#DD0000", width=3),
      "ubiquination"               = list(line="SOLID",          targetArrow="Arrow",    sourceArrow="No Arrow", color="#DD00DD", width=3),
      "repression"                 = list(line="SOLID",          targetArrow="T",        sourceArrow="No Arrow", color="#DD0000", width=3),
      "dissociation"               = list(line="DASH_DOT",       targetArrow="No Arrow", sourceArrow="No Arrow", color="#000000", width=3),
      "indirect effect"            = list(line="EQUAL_DASH",     targetArrow="Arrow",    sourceArrow="No Arrow", color="#000000", width=3),
      "missing interaction"        = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#AA00AA", width=3),
      "expression"                 = list(line="SOLID",          targetArrow="Delta",    sourceArrow="No Arrow", color="#0000AA", width=3),
      "missing interaction"        = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#AA0000", width=3),
      "dimerize"                   = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#00AA00", width=3),
      "phosphorylation"            = list(line="SINEWAVE",       targetArrow="Arrow",    sourceArrow="No Arrow", color="#00AA00", width=3),
      "dephosphorylation"          = list(line="DOT",            targetArrow="Arrow",    sourceArrow="No Arrow", color="#00AA00", width=3),
      "compound"                   = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#000000", width=3),
      "unspecified"                = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#000000", width=3),
      "functional"                 = list(line="DOT",            targetArrow="No Arrow", sourceArrow="No Arrow", color="#000000", width=3)
      )

    if (length(edgeNames(g)) > 0) {
        edgeTypes.in.graph <- unique(as.character(eda(g, "edgeType")))
        unexpected.edgeTypes <- setdiff(edgeTypes.in.graph, names(edge.map))
        if (length(unexpected.edgeTypes) > 0)
           for (type in unexpected.edgeTypes) {
               printf("   unexpected edge type: %s", type);
           } # for type
        stopifnot (all(edgeTypes.in.graph %in% names(edge.map)))
        edge.types <- names(edge.map)

       line.types <- as.character(sapply(edge.map, function(info) info$line))
       setEdgeLineStyleRule(cw, "edgeType", edge.types, line.types)

       target.arrows <- as.character(sapply(edge.map,
                                           function(info) info$targetArrow))
       setEdgeTargetArrowRule(cw, "edgeType", edge.types, target.arrows)

       source.arrows <- as.character(sapply(edge.map,
                                           function(info) info$sourceArrow))
       setEdgeSourceArrowRule(cw, "edgeType", edge.types, source.arrows)

       colors <- as.character(sapply(edge.map, function(info) info$color))
       setEdgeColorRule(cw, "edgeType",  edge.types, colors, mode="lookup")

       setEdgeTargetArrowColorRule(cw, "edgeType",  edge.types, colors)
       setEdgeSourceArrowColorRule(cw, "edgeType",  edge.types, colors)

       widths <- as.integer(sapply(edge.map, function(info) info$width))
       setEdgeLineWidthRule(cw, "edgeType",  edge.types, widths)

       setEdgeTooltipRule(cw, "edgeType")
       } # if at least one edge

    if ("desc" %in% noa.names(cw@graph))
        setNodeTooltipRule(cw, "desc")


    setDefaultNodeFontSize(cw, 60)

    if ("info.node" %in% nodes(cw@graph)) {
        setNodeShapeDirect(cw, "info.node", "rect")
        setNodeSizeDirect(cw, "info.node", 100)
        setNodeFontSizeDirect(cw, "info.node", 120);
        if(!is.null(tumor.name))
           setNodeLabelDirect(cw, "info.node", tumor.name)
        setNodeColorDirect(cw, "info.node", "#FFEEFF")
        setNodeFillOpacityDirect(cw, "info.node", 0)
        setNodeBorderOpacityDirect(cw, "info.node", 0)
        } # info.node

      setTooltipInitialDelay(cw,0)
      setTooltipDismissDelay(cw,100000)

      redraw(cw)
      cw

} # tumorViz
#-------------------------------------------------------------------------------
displayTumor <- function(cw, tumor.name, tbl.mrna, tbl.cnv, tbl.mut)
{
    stopifnot(tumor.name %in% rownames(tbl.mrna))
    z.scores <- as.numeric(tbl.mrna [tumor.name,])
    geneIDs <- colnames(tbl.mrna)
    removers <- which(is.na(z.scores))
    if (length(removers) > 0) {
        z.scores <- z.scores [-removers]
        geneIDs <- geneIDs [-removers]
        }
    setNodeAttributesDirect(cw, "lfc", "numeric", geneIDs, z.scores)
    stopifnot(tumor.name %in% rownames(tbl.cnv))
    copy.number.status <- as.character(tbl.cnv [tumor.name,])
    geneIDs <- colnames(tbl.cnv)
    removers <- which(copy.number.status == "NaN")
    if (length(removers) > 0) {
        copy.number.status <- copy.number.status [-removers]
        geneIDs <- geneIDs [-removers]
        }
    setNodeAttributesDirect(cw, "copyNumber", "char", geneIDs,
                            copy.number.status)

      # determine the mutants
    mutant.gene.list <- extractMutantGenes(tbl.mut, tumor.name)
    geneID.geneSymbol.map <-
      AnnotationDbi::mget(nodes(cw@graph), org.Hs.egSYMBOL, ifnotfound=NA)
    this.tumors.geneID.geneSymbol.map <- geneID.geneSymbol.map
    if (length(mutant.gene.list) > 0) {
         gene.symbols.for.mutants <-
             as.character(geneID.geneSymbol.map [names(mutant.gene.list)])
         mutation.names <- as.character(mutant.gene.list)
         full.mutant.gene.labels <- paste(gene.symbols.for.mutants,
                                         mutation.names)
         this.tumors.geneID.geneSymbol.map [names(mutant.gene.list)] <-
                full.mutant.gene.labels
         }
    setNodeAttributesDirect(cw, "label", "char",
                            names(this.tumors.geneID.geneSymbol.map),
                            as.character(this.tumors.geneID.geneSymbol.map))

      # "score" is a somewhat ad hoc combination of expression, copy number and
      #  mutation, and controls node size(see tumorViz for the rule)

    raw.scores <- as.numeric(tbl.mrna[tumor.name,])
    names(raw.scores) <- colnames(tbl.mrna)

    score <- 4 * abs(raw.scores)

    cnv.genes <- colnames(tbl.cnv)
    homo.deletion.genes <- cnv.genes[which(tbl.cnv [tumor.name, ] == -2)]
    hetero.deletion.genes <- cnv.genes[which(tbl.cnv [tumor.name, ] == -1)]
    gain.genes <- cnv.genes[which(tbl.cnv [tumor.name, ] == 1)]
    amplified.genes <- cnv.genes[which(tbl.cnv [tumor.name, ] == 2)]

    printf("homo.deletion.gene count: %d", length(homo.deletion.genes))
    printf("hetero.deletion.gene count: %d", length(hetero.deletion.genes))
    printf("gain.gene count: %d", length(gain.genes))
    printf("amplified.gene count: %d", length(amplified.genes))

    if (length(homo.deletion.genes) > 0)
        score [homo.deletion.genes] <- score [homo.deletion.genes] + 4

    if (length(hetero.deletion.genes) > 0)
         score [hetero.deletion.genes] <- score [hetero.deletion.genes] + 2

    if (length(gain.genes) > 0)
        score [gain.genes] <- score [gain.genes] + 2

    if (length(amplified.genes) > 0)
         score [amplified.genes] <- score [amplified.genes] + 4

    if (length(mutant.gene.list) > 0)
        score [names(mutant.gene.list)] <-
            score [names(mutant.gene.list)] + 20

    geneIDs <- colnames(tbl.cnv)

    removers <- which(is.na(score))
    if (length(removers) > 0) {
        score <- score [-removers]
        geneIDs <- geneIDs [-removers]
        }
    printf("score: %d", length(score))
    printf("geneIDs: %d", length(geneIDs))
    setNodeAttributesDirect(cw, "score", "numeric", geneIDs, score)

    redraw(cw)
    tumorViz(cw, tumor.name)

    cw

} # displayTumor
#-------------------------------------------------------------------------------
createKeggGraph <- function(tbl, node.name.map, gene.type.map)
{
    gene.ids <- unique(c(tbl$from, tbl$to))
    g <- new("graphNEL", edgemode="directed")
    g <- initEdgeAttribute(g, "edgeType",   "char",    "undefined")
    g <- initEdgeAttribute(g, "source",     "char",    "undefined")
    g <- initNodeAttribute(g, "nodeType",   "char",    "gene")
    g <- initNodeAttribute(g, "label",      "char",    "undefined")
    g <- initNodeAttribute(g, "lfc",        "numeric",  0.0)
    g <- initNodeAttribute(g, "copyNumber", "char",    "0")
    g <- initNodeAttribute(g, "mutation",   "char",    "")
    g <- initNodeAttribute(g, "score",   "numeric",    0.0)
    g <- graph::addNode(gene.ids, g)
    g <- graph::addNode("info.node", g)

    nodeData(g, names(node.name.map), "label") <- as.character(node.name.map)
    nodesWithKnownGeneType <- intersect(names(gene.type.map),
                                        names(node.name.map))
    sub.gene.type.map <- gene.type.map[nodesWithKnownGeneType]
    nodeData(g, names(sub.gene.type.map), "nodeType") <-
              as.character(sub.gene.type.map)

    nodeData(g, "info.node", "nodeType") <- "info"
    nodeData(g, "info.node", "label") <- ""
    g <- graph::addEdge(tbl$from, tbl$to, g)
    edgeData(g, tbl$from, tbl$to, "edgeType") <- tbl$subtype
    edgeData(g, tbl$from, tbl$to, "source") <- tbl$source

    g

} # createKeggGraph
#-------------------------------------------------------------------------------
extractMutantGenes <- function(tbl.mut, tumor.name)
{
    indices <- which(!is.na(tbl.mut [tumor.name,]))
    if (length(indices) == 0)
        return(character(0))
    genes <- colnames(tbl.mut) [indices]
    mutations <- as.character(tbl.mut [tumor.name, indices])
    names(mutations) <- genes

    mutations

} # extractMutantGenes
#-------------------------------------------------------------------------------
dimInactiveNodesAndEdges <- function(cw, tumor.name, genes.of.interest=NULL,
                                    tbl.mrna, tbl.cnv, tbl.mut, mrna.threshold,
                                    opacity)
{
    restoreDefaultOpacities(cw)
    redraw(cw)
    if (is.null(genes.of.interest))
        genes.of.interest <- nodes(cw@graph)
    removers <- grep("info.node", genes.of.interest)
    if (length(removers==1))
        genes.of.interest <- genes.of.interest[-removers]
    active.gene.lists <- determineActiveGenes(tumor.name, genes.of.interest,
                                              tbl.mrna, tbl.cnv, tbl.mut,
                                              mrna.threshold)
    geneIDs.active <- active.gene.lists$all
    geneIDs.extended <- c()
    geneIDs.all <- unique(c(geneIDs.extended, geneIDs.active))
    inactive.geneIDs <- setdiff(nodes(cw@graph), geneIDs.all)
    info.node.index <- grep("info", inactive.geneIDs)
    if (length(info.node.index) > 0)
        inactive.geneIDs <- inactive.geneIDs [-info.node.index]
    setNodeOpacityDirect(cw, inactive.geneIDs, opacity)
    cy2.edges <- getEdgesForNodes(cw@graph, inactive.geneIDs)
    setEdgeOpacityDirect(cw, cy2.edges, opacity)
    redraw(cw)

    inactive.geneIDs

} # dimInactiveNodesAndEdges
#-------------------------------------------------------------------------------
restoreDefaultOpacities <- function(cw)
{
    setEdgeOpacityDirect(cw, cy2.edge.names(cw@graph), 255)
    all.nodes <- nodes(cw@graph)
    info.node.index <- grep("info", all.nodes)
    if (length(info.node.index) > 0)
        all.nodes <- all.nodes[-info.node.index]
    setNodeOpacityDirect(cw, all.nodes, 255)
    redraw(cw)

    cw

} # restoreDefaultOpacities
#-------------------------------------------------------------------------------
determineActiveGenes <- function(tumor.name, genes.of.interest, tbl.mrna,
                                tbl.cnv, tbl.mut, mrna.threshold=0.7)
{
    stopifnot(tumor.name %in% rownames(tbl.mrna))

    mrna.active.indices <-
         which(abs(as.numeric(tbl.mrna [tumor.name,])) > mrna.threshold)
    mrna.genes <- c()
    if (length(mrna.active.indices) > 0)
        mrna.genes <- colnames(tbl.mrna)[mrna.active.indices]

    na.indices <- which(is.na(mrna.genes))
    if (length(na.indices) > 0)
        mrna.genes <- mrna.genes [-na.indices]

    mutant.gene.indices <- which(!is.na(tbl.mut [tumor.name,]))
    mutant.genes <- c()
    if (length(mutant.gene.indices) > 0)
        mutant.genes <- colnames(tbl.mrna)[mutant.gene.indices]

    cnv.indices <- which(abs(tbl.cnv [tumor.name,]) > 0)
    cnv.genes <- c()
    if (length(cnv.indices) > 0)
        cnv.genes <-  colnames(tbl.mrna)[cnv.indices]

    all.active.genes <- unique(c(mrna.genes, mutant.genes, cnv.genes))
    na.indices <- grep("^NA", all.active.genes)

    if (length(na.indices) > 0) {
        all.active.genes <- all.active.genes [-na.indices]
        }

   list(all=intersect(all.active.genes, genes.of.interest),
        mrna=intersect(mrna.genes, genes.of.interest),
        mut=intersect(mutant.genes, genes.of.interest),
        cnv=intersect(cnv.genes, genes.of.interest))

} # determineActiveGenes
#-------------------------------------------------------------------------------
getEdgesForNodes <- function(g, nodeNames)
{
    all.edge.names <- cy2.edge.names(g)
    indices.of.edges.with.nodeNames <- c()
    x <- sapply(nodeNames, function(nodeName)
               grep(nodeName, names(all.edge.names)))
    all.edge.indices <- unlist(x, use.names=FALSE)
        # every edge which index which appears twice must have
        # its a and b nodes both in the nodeNames list
        # find those double hits, get and return those cy2-style edgenames
    edge.tbl <- table(all.edge.indices)
    edge.name.indices <- as.integer(names(which(edge.tbl>=1)))  # == 2

    all.edge.names [edge.name.indices]

} # getEdgesForNodes
#-------------------------------------------------------------------------------


