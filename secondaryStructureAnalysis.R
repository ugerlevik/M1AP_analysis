#########################################################################
## Project: M1AP                                 
## Script purpose: Secondary Structure Analysis
## Date: Dec 19, 2020                                                   
## Author: Umut Gerlevik                                                 
#########################################################################

# Read packages -------------------
library(bio3d)
library(stringi)
library(readr)
library(ggplot2)
library(ggpubr)

# Common variables ----------
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")

# Define functions ----------------
ssa <- function(pdb, traj, start = 1, last = nrow(traj)) {
  require(stringi)
  require(readr)
  require(bio3d)
  
  tmp <- stri_rand_strings(1, 5)
  dir.create(paste0("_outputs/stride_", tmp))
  out <- matrix(ncol = length(pdb[pdb$calpha]))
  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")
  out <- as.data.frame(na.omit(out))
  
  for(i in 1:nrow(traj)) {
    write.pdb(pdb, file = paste0("_outputs/stride_", tmp, "/", tmp, "_", i, ".pdb"), xyz = traj[i,])
    res <- suppressWarnings(system(paste0("stride ", "_outputs/stride_", tmp, "/", tmp, "_", i, ".pdb"), intern = TRUE))
    res <- read_table(as.character(res[grep("ASG", res)]), col_names = FALSE)[, c(2, 3, 4, 6, 7)]
    out <- rbind(out, res$X6)
  }
  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")
  
  unlink(paste0("_outputs/stride_", tmp), recursive = TRUE)
  
  rownames(out) <- start:last
  
  return(as.data.frame(out))
}

ssa_plot <- function(ssa1, ssa2,
                     name1 = "Wild-type", name2,
                     resid1 = as.numeric(gsub('\\D+','', colnames(ssa1)[1])),
                     resid2 = as.numeric(gsub('\\D+','', colnames(ssa1)[length(colnames(ssa1))])),
                     color_number1 = 1, color_number2 = 2) {
  palet <- ggsci::pal_jama()
  colorPalet <- palet(7)
  
  resids <- resid1:resid2
  
  out <- matrix(ncol = 7, nrow = length(resids))
  # H: AlphaHelix, E: ExtendedConformation(BetaSheet), B or b: Bridge,
  # T: Turn, C or " ": Coil, G: Helix310, I: PiHelix
  colnames(out) <- c("AlphaHelix", "BetaSheet", "Bridge", "Turn",
                     "Coil", "Helix310", "PiHelix")
  
  out1 <- as.data.frame(na.omit(out))
  for(i in resids) {
    perc <- (table(ssa1[i]) / sum(table(ssa1[i])))*100
    
    out1[i, "AlphaHelix"] <- perc["H"]
    out1[i, "BetaSheet"] <- perc["E"]
    out1[i, "Bridge"] <- sum(perc[c("B","b")])
    out1[i, "Turn"] <- perc["T"]
    out1[i, "Coil"] <- perc["C"]
    out1[i, "Helix310"] <- perc["G"]
    out1[i, "PiHelix"] <- perc["I"]
  }
  out2 <- as.data.frame(na.omit(out))
  for(i in resids) {
    perc <- (table(ssa2[i]) / sum(table(ssa2[i])))*100
    
    out2[i, "AlphaHelix"] <- perc["H"]
    out2[i, "BetaSheet"] <- perc["E"]
    out2[i, "Bridge"] <- sum(perc[c("B","b")])
    out2[i, "Turn"] <- perc["T"]
    out2[i, "Coil"] <- perc["C"]
    out2[i, "Helix310"] <- perc["G"]
    out2[i, "PiHelix"] <- perc["I"]
  } 
  
  out1 <- out1[resid1:resid2, ]
  out1$resid <- resids
  out1 <- reshape2::melt(out1[1:8], id.var = "resid")
  out1$mut <- name1
  
  out2 <- out2[resid1:resid2, ]
  out2$resid <- resids
  out2 <- reshape2::melt(out2[1:8], id.var = "resid")
  out2$mut <- name2
  
  out <- rbind(out1, out2)
  
  require(ggplot2)
  require(ggpubr)
  p1 <- ggplot(out[out$variable == "AlphaHelix",], aes(x = resid, y = value, fill = mut,)) 
  p1 <- p1 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p1 <- p1 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p1 <- p1 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p1 <- p1 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p1 <- p1 + labs(x = "Residue Index", y = "Alpha Helix (%)")
  p1 <- p1 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p1 <- p1 + theme_minimal()
  p1 <- p1 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p2 <- ggplot(out[out$variable == "BetaSheet",], aes(x = resid, y = value, fill = mut,)) 
  p2 <- p2 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p2 <- p2 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p2 <- p2 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p2 <- p2 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p2 <- p2 + labs(x = "Residue Index", y = "Beta Sheet (%)")
  p2 <- p2 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p2 <- p2 + theme_minimal()
  p2 <- p2 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p3 <- ggplot(out[out$variable == "Bridge",], aes(x = resid, y = value, fill = mut,)) 
  p3 <- p3 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p3 <- p3 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p3 <- p3 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p3 <- p3 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p3 <- p3 + labs(x = "Residue Index", y = "Beta Bridge (%)")
  p3 <- p3 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p3 <- p3 + theme_minimal()
  p3 <- p3 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p4 <- ggplot(out[out$variable == "Turn",], aes(x = resid, y = value, fill = mut,)) 
  p4 <- p4 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p4 <- p4 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p4 <- p4 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p4 <- p4 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p4 <- p4 + labs(x = "Residue Index", y = "Turn (%)")
  p4 <- p4 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p4 <- p4 + theme_minimal()
  p4 <- p4 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p5 <- ggplot(out[out$variable == "Coil",], aes(x = resid, y = value, fill = mut,)) 
  p5 <- p5 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p5 <- p5 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p5 <- p5 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p5 <- p5 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p5 <- p5 + labs(x = "Residue Index", y = "Coil (%)")
  p5 <- p5 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p5 <- p5 + theme_minimal()
  p5 <- p5 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p6 <- ggplot(out[out$variable == "Helix310",], aes(x = resid, y = value, fill = mut,)) 
  p6 <- p6 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p6 <- p6 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p6 <- p6 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p6 <- p6 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p6 <- p6 + labs(x = "Residue Index", y = "3(10) Helix (%)")
  p6 <- p6 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p6 <- p6 + theme_minimal()
  p6 <- p6 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p7 <- ggplot(out[out$variable == "PiHelix",], aes(x = resid, y = value, fill = mut,)) 
  p7 <- p7 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p7 <- p7 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p7 <- p7 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p7 <- p7 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p7 <- p7 + labs(x = "Residue Index", y = "Pi Helix (%)")
  p7 <- p7 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p7 <- p7 + theme_minimal()
  p7 <- p7 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_text(angle = 90, size = 10),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
  
  p8 <- ggpubr::as_ggplot(ggpubr::get_legend(p1 + theme(legend.position = "top")))
  
  
  fig <- ggarrange(p8, p1, p2, p3, p4, p5, p6, p7, ncol = 1,
                   heights = c(1, rep(5, 7)))
  
  return(annotate_figure(fig, bottom = text_grob("Residue Index", face = "bold", size = 18)))
}


# Calculate Secondary Structures --------------
# pdb_WT <- read.pdb("3_WT/3_production/seokWT_onlyProt.pdb")
# dcd_WT_r1 <- read.dcd(trjfile = "3_WT/3_production/seokWT_250ns_onlyProt_each25k.dcd")
# dcd_WT_r2 <- read.dcd(trjfile = "3_WT/4_production_repeat2/seokWT_250ns_r2_onlyProt_each25k.dcd")
# ssa_WT <- ssa(pdb_WT, rbind(dcd_WT_r1, dcd_WT_r2))
# saveRDS(ssa_WT, "_outputs/ssa/ssa_WT.RDS")
# rm(pdb_WT, dcd_WT_r1[seq(1, nrow(dcd_WT_r1), 2),], dcd_WT_r2[seq(1, nrow(dcd_WT_r2), 2),])

# pdb_S50P <- read.pdb("7_S50P/3_production/seokS50P_onlyProt.pdb")
# dcd_S50P_r1 <- read.dcd(trjfile = "7_S50P/3_production/seokS50P_250ns_onlyProt_each25k.dcd")
# dcd_S50P_r2 <- read.dcd(trjfile = "7_S50P/4_production_repeat2/seokS50P_250ns_r2_onlyProt_each25k.dcd")
# ssa_S50P <- ssa(pdb_S50P,
#                 rbind(dcd_S50P_r1[seq(1, nrow(dcd_S50P_r1), 2),],dcd_S50P_r2[seq(1, nrow(dcd_S50P_r2), 2),]))
# saveRDS(ssa_S50P, "_outputs/ssa/ssa_S50P.RDS")
# rm(pdb_S50P, dcd_S50P_r1, dcd_S50P_r2)

# pdb_R266Q <- read.pdb("8_R266Q/3_production/seokR266Q_onlyProt.pdb")
# dcd_R266Q_r1 <- read.dcd(trjfile = "8_R266Q/3_production/seokR266Q_250ns_onlyProt_each25k.dcd")
# dcd_R266Q_r2 <- read.dcd(trjfile = "8_R266Q/4_production_repeat2/seokR266Q_250ns_r2_onlyProt_each25k.dcd")
# ssa_R266Q <- ssa(pdb_R266Q,
#                  rbind(dcd_R266Q_r1[seq(1, nrow(dcd_R266Q_r1), 2),], dcd_R266Q_r2[seq(1, nrow(dcd_R266Q_r2), 2),]))
# saveRDS(ssa_R266Q, "_outputs/ssa/ssa_R266Q.RDS")
# rm(pdb_R266Q, dcd_R266Q_r1, dcd_R266Q_r2)

# pdb_G317R <- read.pdb("5_G317R/3_production/seokG317R_onlyProt.pdb")
# dcd_G317R_r1 <- read.dcd(trjfile = "5_G317R/3_production/seokG317R_250ns_onlyProt_each25k.dcd")
# dcd_G317R_r2 <- read.dcd(trjfile = "5_G317R/4_production_repeat2/seokG317R_250ns_r2_onlyProt_each25k.dcd")
# ssa_G317R <- ssa(pdb_G317R, rbind(dcd_G317R_r1[seq(1, nrow(dcd_G317R_r1), 2),], dcd_G317R_r2[seq(1, nrow(dcd_G317R_r2), 2),]))
# saveRDS(ssa_G317R, "_outputs/ssa/ssa_G317R.RDS")
# rm(pdb_G317R, dcd_G317R_r1, dcd_G317R_r2)

# pdb_P389L <- read.pdb("4_P389L/3_production/seokP389L_onlyProt.pdb")
# dcd_P389L_r1 <- read.dcd(trjfile = "4_P389L/3_production/seokP389L_250ns_onlyProt_each25k.dcd")
# dcd_P389L_r2 <- read.dcd(trjfile = "4_P389L/4_production_repeat2/seokP389L_250ns_r2_onlyProt_each25k.dcd")
# ssa_P389L <- ssa(pdb_P389L, rbind(dcd_P389L_r1[seq(1, nrow(dcd_P389L_r1), 2),], dcd_P389L_r2[seq(1, nrow(dcd_P389L_r2), 2),]))
# saveRDS(ssa_P389L, "_outputs/ssa/ssa_P389L.RDS")
# rm(pdb_P389L, dcd_P389L_r1, dcd_P389L_r2)

# pdb_L430P <- read.pdb("6_L430P/3_production/seokL430P_onlyProt.pdb")
# dcd_L430P_r1 <- read.dcd(trjfile = "6_L430P/3_production/seokL430P_250ns_onlyProt_each25k.dcd")
# dcd_L430P_r2 <- read.dcd(trjfile = "6_L430P/4_production_repeat2/seokL430P_250ns_r2_onlyProt_each25k.dcd")
# ssa_L430P <- ssa(pdb_L430P, rbind(dcd_L430P_r1[seq(1, nrow(dcd_L430P_r1), 2),], dcd_L430P_r2[seq(1, nrow(dcd_L430P_r2), 2),]))
# saveRDS(ssa_L430P, "_outputs/ssa/ssa_L430P.RDS")
# rm(pdb_L430P, dcd_L430P_r1, dcd_L430P_r2)

# Visualize ------------------------------
ssa_WT <- readRDS("_outputs/ssa/ssa_WT.RDS")
ssa_S50P <- readRDS("_outputs/ssa/ssa_S50P.RDS")
ssa_R266Q <- readRDS("_outputs/ssa/ssa_R266Q.RDS")
ssa_G317R <- readRDS("_outputs/ssa/ssa_G317R.RDS")
ssa_P389L <- readRDS("_outputs/ssa/ssa_P389L.RDS")
ssa_L430P <- readRDS("_outputs/ssa/ssa_L430P.RDS")

p_S50P <- ssa_plot(ssa_WT, ssa_S50P, name2 = mutationlist[1], color_number2 = 2)
p_R266Q <- ssa_plot(ssa_WT, ssa_R266Q, name2 = mutationlist[2], color_number2 = 3)
p_G317R <- ssa_plot(ssa_WT, ssa_G317R, name2 = mutationlist[3], color_number2 = 4)
p_P389L <- ssa_plot(ssa_WT, ssa_P389L, name2 = mutationlist[4], color_number2 = 5)
p_L430P <- ssa_plot(ssa_WT, ssa_L430P, name2 = mutationlist[5], color_number2 = 6)

pdf("_outputs/secondaryStructure.pdf", width = 20, height = 20)
p_S50P
p_R266Q
p_G317R
p_P389L
p_L430P
dev.off()

# focused:
p_focused_S50P <- ssa_plot(ssa_WT, ssa_S50P, name2 = mutationlist[1], color_number2 = 2, resid1 = 45, resid2 = 55)
p_focused_R266Q <- ssa_plot(ssa_WT, ssa_R266Q, name2 = mutationlist[2], color_number2 = 3, resid1 = 261, resid2 = 271)
p_focused_G317R <- ssa_plot(ssa_WT, ssa_G317R, name2 = mutationlist[3], color_number2 = 4, resid1 = 312, resid2 = 322)
p_focused_P389L <- ssa_plot(ssa_WT, ssa_P389L, name2 = mutationlist[4], color_number2 = 5, resid1 = 384, resid2 = 394)
p_focused_L430P <- ssa_plot(ssa_WT, ssa_L430P, name2 = mutationlist[5], color_number2 = 6, resid1 = 425, resid2 = 435)
empty <- ggplot() + theme_void()

fig <- ggarrange(ggarrange(p_focused_S50P, p_focused_R266Q, p_focused_G317R,
                           ncol = 3, nrow = 1, align = "h", widths = c(1, 1, 1),
                           labels = c("A", "B", "C"), font.label = list(size = 28, face = "bold")),
                 ggarrange(empty, p_focused_P389L, p_focused_L430P, empty,
                           ncol = 4, nrow = 1, align = "h", widths = c(.5, 1, 1, .5),
                           labels = c("", "D", "E", ""), font.label = list(size = 28, face = "bold")),
                 nrow = 2, ncol = 1, widths = c(1.5, 1))

pdf("_outputs/secondaryStructure_focused.pdf", width = 20, height = 45)
fig
dev.off()
