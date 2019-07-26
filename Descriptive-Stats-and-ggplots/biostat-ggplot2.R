## ----exDF----------------------------------------------------------------
people <- data.frame(weight = c(80, 49, 62, 57),
                     height = c(1.82, 1.58, 1.71, 1.63),
                     gender = c("M", "F", "F", "F"))
people

## ----grammarPlot, fig.width=15, echo=FALSE-------------------------------
library(ggplot2)
library(gridExtra)
p1 <- ggplot(people, aes(y = weight, x = height, colour = gender)) +
  geom_point(size = 8) + theme_void() + guides(colour = FALSE)
p2 <- ggplot(people, aes(y = weight, x = height, colour = gender)) +
  geom_point(shape = NA) + guides(colour = FALSE) + theme_bw() +
  theme(text = element_text(size = 20))
p3 <- ggplot(people, aes(y = weight, x = height, colour = gender)) +
  geom_point(size = 8) + guides(colour = FALSE) + theme_bw() +
  theme(text = element_text(size = 20))
grid.arrange(p1, p2, p3, ncol = 3)

## ---- eval=FALSE---------------------------------------------------------
## library(ggplot2)
## ggplot(data = <DATA>,
##        aes(x = <X AXIS VARIABLE>,
##            y = <Y AXIS VARIABLE>, ... ), ...) +
##
##   geom_<TYPE>(aes(size = <SIZE VARIABLE>, ...),
##                    data = <DATA>,
##                    stat = <FUNCTION>,
##                    position = <POSITION>,
##                    color = <"COLOR">, ...) +
##
##   scale_<AESTHETIC>_<TYPE>(name = <NAME>,
##                    breaks = <WHERE>,
##                    labels = <LABELS>, ... ) +
##
##   theme(...) +
##   facet_<TYPE>(<FORMULA>)

## ----diamondDF-----------------------------------------------------------
data(diamonds)
dim(diamonds)
head(diamonds)

## ----subsetDiam, echo=FALSE----------------------------------------------
set.seed(25091309)
sample1000 <- sample(1:nrow(diamonds), 1000, replace = FALSE)
diamonds <- diamonds[sample1000, ]

## ----baseHist------------------------------------------------------------
hist(diamonds$price, main = "", xlab = "Price", breaks = 50)

## ----ggplotHist----------------------------------------------------------
ggplot(diamonds, aes(x = price)) + geom_histogram(bins = 50)

## ----ggplotHist2---------------------------------------------------------
ggplot(diamonds, aes(x = price, fill = color)) +
  geom_histogram(binwidth = 1000) + facet_wrap( ~ cut) + theme_bw()

## ----breaking, fig.width=12, fig.height=4--------------------------------
## Create ggplot object, populate it with data
ggplot(diamonds, aes(x = carat, y = price, colour = cut)) +

## Add layer(s)
  geom_point(alpha = 0.3) +
  geom_smooth() +

## Scales for dimensions, color palettes
	scale_y_log10() +

## Condition on variables
	facet_grid(~ cut) +

## More options
	ggtitle("First example") + theme_bw()

## ----breaking2, eval=FALSE-----------------------------------------------
## ## Create ggplot object
## MyPlot <- ggplot(diamonds)
## class(MyPlot)
## summary(MyPlot); MyPlot
##
## ## Add aesthetics
## MyPlot <- MyPlot + aes(x = carat, y = price, colour = cut)
## summary(MyPlot)
## MyPlot
##
## ## Add layer(s)
## MyPlot <- MyPlot + geom_point(alpha=0.3)
## summary(MyPlot) ; MyPlot
##
## MyPlot <- MyPlot + geom_smooth()
## summary(MyPlot) ; MyPlot
##
## ## Scales for dimensions
## MyPlot + scale_y_log10()
##
## ## Condition on variables
## MyPlot + facet_grid(~ cut)
##
## ## More options
## MyPlot + ggtitle("First example") + theme_bw()

## ----helpAES, eval=FALSE-------------------------------------------------
## ?geom_point
## ...
## Aesthetics
## The following aesthetics can be used with geom_point. Aesthetics are mapped
## to variables in the data with the aes function: geom_point(aes(x = var))
## x: x position (required)
## y: y position (required)
## shape: shape of point
## colour: border colour
## size: size
## fill: internal colour
## alpha: transparency

## ----equivAES, eval=FALSE------------------------------------------------
## ggplot(diamonds, aes(x = carat, y = price, color = cut)) + geom_point()
##
## ggplot(diamonds) +  geom_point(aes(x = carat, y = price, color = cut))
##
## ggplot(diamonds, aes(x = carat, y = price)) + geom_point(aes(color = cut))

## ----equivAESrun, echo=FALSE, fig.height=5-------------------------------
ggplot(diamonds, aes(x = carat, y = price, color = cut)) + geom_point()

## ----AESvariable, eval=FALSE---------------------------------------------
## ggplot(diamonds, aes(x = carat, y = price, color = clarity)) + geom_point()

## ----AESfixed, fig.height=4----------------------------------------------
ggplot(diamonds, aes(x = carat, y = price)) + geom_point(color = "red")

## ----ex1, echo=FALSE-----------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price, shape = cut)) + geom_point()

## ----ex1b, echo=FALSE----------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price, shape = cut)) + geom_point()

## ----ex2, echo=FALSE-----------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price, shape = cut)) +
  geom_point(alpha = 0.2)

## ----ex3, echo=FALSE-----------------------------------------------------
ggplot(diamonds, aes(x = carat, fill = cut)) + geom_histogram(binwidth = 0.2)

## ----helpGeom, eval=FALSE------------------------------------------------
## help.search("geom_", package = "ggplot2")

## ---- eval=FALSE---------------------------------------------------------
## geom_abline       geom_jitter
## geom_area		      geom_line
## geom_bar		      geom_linerange
## geom_bin2d		    geom_path
## geom_blank		    geom_point
## geom_boxplot	    geom_pointrange
## geom_contour	    geom_polygon
## geom_crossbar	    geom_quantile
## geom_density	    geom_rect
## geom_density2d	  geom_ribbon
## geom_errorbar	    geom_rug
## geom_errorbarh	  geom_segment
## geom_freqpoly		  geom_smooth
## geom_hex		      geom_step
## geom_histogram	  geom_text
## geom_hline		    geom_tile
## geom_vline        ...

## ----hist1---------------------------------------------------------------
p <- ggplot(diamonds)
## Overall histogram
p + geom_histogram(aes(x = price))

## ----hist2---------------------------------------------------------------
## Composition of each bin
p + geom_histogram(aes(x = price, fill = cut))

## ----hist3---------------------------------------------------------------
## Relative proportions
p + geom_histogram(aes(x = price, fill = cut), position = "fill")

## ----dens----------------------------------------------------------------
p + geom_density(aes(x = price, fill = cut), alpha = 0.5)

## ----box1----------------------------------------------------------------
p + geom_boxplot(aes(x = cut, y = price), notch = TRUE)

## ----box2----------------------------------------------------------------
p + geom_boxplot(aes(x = cut, y = price, fill = color))

## ----rect, fig.height=4--------------------------------------------------
p + geom_tile(aes(x = as.numeric(cut), y = as.numeric(color), fill = depth))

## ----rect2---------------------------------------------------------------
p + geom_raster(aes(x = as.numeric(cut), y = as.numeric(color), fill = depth))

## ----multi1--------------------------------------------------------------
ggplot(diamonds, aes(x = color, y = price, fill = color)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(aes(fill = color), alpha = 0.1, shape = 21)

## ----multi2--------------------------------------------------------------
ggplot(diamonds, aes(x = color, y = price, fill = color)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(aes(fill = color), alpha = 0.1, , shape = 21,
             position = position_jitter(w = .3))

## ----multi3, fig.height=4------------------------------------------------
ggplot(diamonds, aes(x = reorder(color, price), y = price, fill = color)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(aes(fill = color), alpha = 0.1, , shape = 21,
             position = position_jitter(w = .4))

## ----multi3b-------------------------------------------------------------
ggplot(diamonds, aes(x = reorder(color, price), y = price, fill = color)) +
  geom_violin() +
  geom_point(aes(fill = color), alpha = 0.5, , shape = 21,
             position = position_jitter(w = .4))

## ----multi4--------------------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price, color = cut)) +
  geom_point(shape = 21) + geom_smooth() + geom_rug()

## ----ex4, echo=FALSE, fig.width=5, fig.height=5--------------------------
ggplot(diamonds, aes(x = carat, y = price, colour = depth)) +
  geom_point(alpha = 0.5) + geom_smooth(colour = "darkred")

## ----ex4b, echo=FALSE, fig.width=5, fig.height=5-------------------------
ggplot(diamonds, aes(x = carat, y = price, colour = depth)) +
  geom_point(alpha = 0.5) + geom_vline(colour = "green", xintercept = 2)

## ----ex5, echo=FALSE, fig.width=5, fig.height=5--------------------------
ggplot(diamonds, aes(x = cut, fill = color)) + geom_bar()

## ----ex7, echo=FALSE, fig.width=5, fig.height=5--------------------------
ggplot(diamonds, aes(x = cut, fill = color)) +
  geom_bar(position = position_dodge())

## ----ex6, echo=FALSE, fig.width=5, fig.height=5--------------------------
ggplot(diamonds, aes(x = cut, fill = color)) +
  geom_bar(position = position_fill())

## ----theme1--------------------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price)) + geom_point() +
  geom_smooth(aes(colour = cut)) + theme_bw()

## ----labels--------------------------------------------------------------
ggplot(diamonds, aes(x = carat, y = price)) + geom_point() +
  geom_smooth(aes(colour = cut)) + theme_bw() + ylab("price (in USD)") +
  ggtitle("My beautiful plot")

## ----legend1-------------------------------------------------------------
ggplot(diamonds, aes(x = cut, y = price, fill = clarity)) +
  geom_boxplot() + scale_fill_discrete(name = "Clarity of diamond",
                                       labels = paste("C", 1:8))

## ----legend2-------------------------------------------------------------
library(RColorBrewer)
ggplot(diamonds, aes(x = cut, y = price, fill = clarity)) +
  geom_boxplot() + scale_fill_manual(name = "Clarity of\n diamond",
                                     labels = paste("C", 1:8),
                                     values = brewer.pal(8, "Set2"))

## ----legend3-------------------------------------------------------------
ggplot(diamonds, aes(x = z, y = carat, colour = price)) +
  geom_point() + scale_y_log10() +
  scale_colour_gradient(low = "grey", high = "pink")

## ----legend4-------------------------------------------------------------
ggplot(diamonds, aes(x = z, y = carat, colour = price)) +
  geom_point() + scale_y_log10() + xlim(0, 10) +
  scale_colour_gradient(low = "grey", high = "pink")

## ----facet1--------------------------------------------------------------
ggplot(diamonds, aes(x = price)) + geom_histogram() + facet_wrap(~ cut)

## ----facet2--------------------------------------------------------------
ggplot(diamonds, aes(x = price)) + geom_histogram(fill = "red") +
  facet_grid(clarity ~ cut) + theme_dark()

## ----custom--------------------------------------------------------------
ggplot(diamonds, aes(x = cut, y = price, fill = clarity)) +
  geom_boxplot() +
  theme(legend.text = element_text(size = 5, colour = "red"),
        legend.position = "top",
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, face = "bold"))

## ----paletteBrewer, echo=FALSE-------------------------------------------
ggplot(diamonds, aes(x = x, y = depth, colour = clarity)) +
  geom_point() + scale_color_brewer() + xlab("dimension") +
  ggtitle("How is clarity affected by dimensions?")

## ----paletteBrewer2, echo=FALSE------------------------------------------
ggplot(diamonds, aes(x = x, y = depth, colour = clarity)) +
  geom_point() + scale_color_brewer(type = "qual") + xlab("dimension") +
  ggtitle("How is clarity affected by dimensions?")

## ----bewerPanel, echo=FALSE----------------------------------------------
ggplot(diamonds, aes(x = cut, y = carat, fill = cut)) +
  geom_boxplot() + scale_fill_brewer(name = "cut quality", palette = 7) +
  facet_wrap(~ color) + theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## ----bewerPanel2, echo=FALSE---------------------------------------------
ggplot(diamonds, aes(x = cut, y = carat, fill = cut)) +
  geom_boxplot() + scale_fill_brewer(type = "div", palette = 3) +
  xlab("cut quality") + facet_grid(color ~ clarity) + theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 5, angle = 45))

## ----gridextra, fig.width=18, fig.height=6-------------------------------
library(gridExtra)
p1 <- ggplot(diamonds) + geom_point(aes(x = carat, y = price, color = cut))
p2 <- ggplot(diamonds) + geom_density(aes(x = price, fill = cut), alpha = 0.5)
grid.arrange(p1, p2, ncol = 2)

## ----ggsave, eval=FALSE--------------------------------------------------
## p <- ggplot(...) + ...
## ggsave("...", plot = p, width = 4, height = 4)

## ----ggnetwork, fig.width=6, fig.height=6--------------------------------
library(ggnetwork); library(network)
data(emon)
ggplot(ggnetwork(emon[[1]], layout = "kamadakawai", arrow.gap = 0.025),
  aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(color = Frequency), curvature = 0.1,
  arrow = arrow(length = unit(10, "pt"), type = "open")) +
  geom_nodes(aes(size = Formalization)) +
  scale_color_gradient(low = "grey50", high = "tomato") +
  scale_size_area(breaks = 1:3) + theme_blank()

