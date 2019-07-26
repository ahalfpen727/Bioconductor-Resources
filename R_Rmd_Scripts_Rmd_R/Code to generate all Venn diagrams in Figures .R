library(VennDiagram);

# Figure 1A
venn.diagram(
	x = list(
		Label = 1:100
		),
	filename = "1A-single_Venn.tiff",
	col = "black",
	lwd = 9,
	fontface = "bold",
	fill = "grey",
	alpha = 0.75,
	cex = 4,
	cat.cex = 3,
	cat.fontface = "bold",
	);

# Figure 1B
venn.diagram(
	x = list(
		X = 1:150,
		Y = 121:180
		),
	filename = "1B-double_Venn.tiff",
	lwd = 4,
	fill = c("cornflowerblue", "darkorchid1"),
	alpha = 0.75,
	label.col = "white",
	cex = 4,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("cornflowerblue", "darkorchid1"),
	cat.cex = 3,
	cat.fontfamily = "serif",
	cat.fontface = "bold",
	cat.dist = c(0.03, 0.03),
	cat.pos = c(-20, 14)
	);

# Figure 1C
venn.diagram(
	x = list(
		R = c(1:70, 71:110, 111:120, 121:140),
		B = c(141:200, 71:110, 111:120, 201:230),
		G = c(231:280, 111:120, 121:140, 201:230)
		),
	filename = "1C-triple_Venn.tiff",
	col = "transparent",
	fill = c("red", "blue", "green"),
	alpha = 0.5,
	label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col = c("darkred", "darkblue", "darkgreen"),
	cat.cex = 2.5,
	cat.fontfamily = "serif",
	cat.dist = c(0.06, 0.06, 0.03),
	cat.pos = 0
	);

# Figure 1D
venn.diagram(
	x = list(
		I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205, 206:220),
		IV = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175, 176:180),
		II = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375, 376:405),
		III = c(406:475, 286:335, 106:140, 141:160, 166:175, 181:205, 336:375, 476:530)
		),
	filename = "1D-quadruple_Venn.tiff",
	col = "black",
	lty = "dotted",
	lwd = 4,
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
	alpha = 0.50,
	label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
	cat.cex = 2.5,
	cat.fontfamily = "serif"
	);

# Figure 2-1
venn.diagram(
	x = list(
		A = 1:105,
		B = 101:115
		),
	filename = "2-1_special_case_ext-text.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = c(-20, 20),
	ext.line.lty = "dotted",
	ext.line.lwd = 2,
	ext.pos = 12,
	ext.dist = -0.12,
	ext.l = 0.85
	);

# Figure 2-2
venn.diagram(
	x = list(
		A = 1:100,
		B = 1:10
		),
	filename = "2-2_special_case_pairwise-inclusion.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = 0
	);

# Figure 2-3
venn.diagram(
	x = list(
		A = 1:150,
		B = 151:250
		),
	filename = "2-3_special_case_pairwise-exclusion.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = c(0, 0),
	cat.dist = 0.05
	);

# Figure 2-4
venn.diagram(
	x = list(
		A = c(1:50, 101:140, 141:160, 161:170),
		B = c(171:230, 101:140, 161:170, 291:320),
		C = c(141:160, 161:170, 291:320)
		),
	sp.cases = TRUE,
	filename = "2-4_triple_special_case-001.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.dist = c(0.05, 0.05, -0.1)
	);

# Figure 2-5
venn.diagram(
	x = list(
		A = c(1:100),
		B = c(61:70, 71:100),
		C = c(41:60, 61:70)
		),
	sp.cases = TRUE,
	filename = "2-5_triple_special_case-012AA.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = c(-25, 0, 30),
	cat.dist = c(0.05, 0.05, 0.02)
	);

# Figure 2-6
venn.diagram(
	x = list(
		A = c(1:90),
		B = c(1:25),
		C = c(1:5)
		),
	sp.cases = TRUE,
	filename = "2-6_triple_special_case-022AAAO.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = 0,
	cat.dist = c(0.65, 0.65, 0.15)
	);

# Figure 2-7
venn.diagram(
	x = list(
		A = c(1:20),
		B = c(21:80),
		C = c(81:210)
		),
	sp.cases = TRUE,
	filename = "2-7_triple_special_case-100.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.dist = 1.5
	);

# Figure 2-8
venn.diagram(
	x = list(
		A = c(1:80),
		B = c(41:150),
		C = c(71:100)
		),
	sp.cases = TRUE,
	filename = "2-8_triple_special_case-011A.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.dist = c(0.07, 0.07, 0.02),
	cat.pos = c(-20, 20, 20)
	);

# Figure 2-9
venn.diagram(
	x = list(
		A = c(1:10),
		B = c(11:90),
		C = c(81:90)
		),
	sp.cases = TRUE,
	filename = "2-9_triple_special_case-121AO.tiff",
	cex = 2.5,
	cat.cex = 2.5,
	cat.pos = 0,
	cat.dist = c(1, 1, 0.3),
	reverse = TRUE
	);
