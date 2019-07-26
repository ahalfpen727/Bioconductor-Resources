## ---- echo = FALSE--------------------------------------------------------------------------------
options(width = 100)

## -------------------------------------------------------------------------------------------------
library(iCOBRA)
data(cobradata_example)
class(cobradata_example)
cobradata_example

## -------------------------------------------------------------------------------------------------
cobradata_example <- calculate_adjp(cobradata_example)

## -------------------------------------------------------------------------------------------------
cobraperf <- calculate_performance(cobradata_example, binary_truth = "status", 
                                  cont_truth = "logFC", splv = "none",
                                  maxsplit = 4)
slotNames(cobraperf)

## -------------------------------------------------------------------------------------------------
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", 
                                  facetted = TRUE)

## ---- fig.width = 7, fig.height = 4, fig.cap = "", warning = FALSE--------------------------------
plot_tpr(cobraplot)

## ---- fig.width = 7, fig.height = 4, fig.cap = "", warning = FALSE--------------------------------
plot_fdrtprcurve(cobraplot)

## ---- fig.width = 7, fig.height = 4, fig.cap = "", warning = FALSE--------------------------------
plot_overlap(cobraplot)

## ---- fig.width = 7, fig.height = 5, warning = FALSE----------------------------------------------
cobraperf <- calculate_performance(cobradata_example, binary_truth = "status", 
                                  cont_truth = "status", splv = "expr_cat")
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", 
                                  facetted = TRUE)
plot_tpr(cobraplot)

## ---- fig.width = 7, fig.height = 5, warning = FALSE----------------------------------------------
library(ggplot2)
pp <- plot_tpr(cobraplot, stripsize = 7.5, pointsize = 3)
pp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                      hjust = 1, size = 10),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 10),
           axis.title.y = element_text(size = 10))

## ---- fig.width = 7, fig.height = 4, warning = FALSE----------------------------------------------
pp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                      hjust = 1, size = 10),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 10),
           axis.title.y = element_text(size = 10),
           legend.position = "bottom") +
  facet_wrap(~splitval, nrow = 1)

## ---- fig.width = 7, fig.height = 5, warning = FALSE----------------------------------------------
plot_overlap(cobraplot)
plot_overlap(cobraplot, cex = c(1, 0.7, 0.7))

## ---- fig.width = 7, fig.height = 5, warning = FALSE----------------------------------------------
cobraplot <- prepare_data_for_plot(cobraperf, 
                                  colorscheme = c("blue", "green", "pink"), 
                                  facetted = TRUE)
pp <- plot_tpr(cobraplot, stripsize = 7.5, pointsize = 3)
pp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                      hjust = 1, size = 10),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 10),
           axis.title.y = element_text(size = 10))

## ---- eval = FALSE--------------------------------------------------------------------------------
#  COBRAapp(cobradata_example)
#  COBRAapp()

## -------------------------------------------------------------------------------------------------
COBRAData_to_text(cobradata = cobradata_example, 
                  truth_file = "cobradata_truth.txt",
                  result_files = "cobradata_results.txt", 
                  feature_id = "feature")

## -------------------------------------------------------------------------------------------------
cobra <- COBRAData_from_text(truth_file = "cobradata_truth.txt", 
                             result_files = "cobradata_results.txt", 
                             feature_id = "feature")
cobra

