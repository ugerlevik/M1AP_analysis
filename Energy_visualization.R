##################################################
## Project: M1AP mutations
## Script purpose: Visualization of Energy
## Date: Dec 30, 2021
## Author: Umut Gerlevik
##################################################

library(ggplot2)
library(gg.gap)
library(ggpubr)

# columns <- unlist(as.list(data.table::fread(text = readLines("3_WT/1_eq1/seokWT_eq1.log")[grep("^ETITLE", readLines("3_WT/1_eq1/seokWT_eq1.log"))])[1]))[2:length(columns)]
# names(columns) <- NULL
# dput(columns)
columns <- c("TS", "BOND", "ANGLE", "DIHED", "IMPRP", "ELECT", "VDW", "BOUNDARY", 
             "MISC", "KINETIC", "TOTAL", "TEMP", "POTENTIAL", "TOTAL3", "TEMPAVG", 
             "PRESSURE", "GPRESSURE", "VOLUME", "PRESSAVG", "GPRESSAVG")

# Wild-type --------------------------------------
sel <- grep("^ENERGY", readLines("3_WT/1_eq1/seokWT_eq1.log"))
WT_eq1 <- data.table::fread(text = readLines("3_WT/1_eq1/seokWT_eq1.log")[sel])[, -1]
colnames(WT_eq1) <- columns

sel <- grep("^ENERGY", readLines("3_WT/2_eq2/seokWT_eq2.log"))
WT_eq2 <- data.table::fread(text = readLines("3_WT/2_eq2/seokWT_eq2.log")[sel])[, -1]
colnames(WT_eq2) <- columns

sel <- grep("^ENERGY", readLines("3_WT/3_production/seokWT_50ns.log"))
WT_50ns <- data.table::fread(text = readLines("3_WT/3_production/seokWT_50ns.log")[sel])[, -1]
colnames(WT_50ns) <- columns
WT_50ns$TS <- WT_50ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_50ns <- WT_50ns[WT_50ns$TS %in% seq(WT_50ns$TS[1], WT_50ns$TS[length(WT_50ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("3_WT/3_production/seokWT_50to100ns.log"))
WT_50to100ns <- data.table::fread(text = readLines("3_WT/3_production/seokWT_50to100ns.log")[sel])[, -1]
colnames(WT_50to100ns) <- columns
WT_50to100ns$TS <- WT_50to100ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_50to100ns <- WT_50to100ns[WT_50to100ns$TS %in% seq(WT_50to100ns$TS[1], WT_50to100ns$TS[length(WT_50to100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("3_WT/3_production/slurm-3548241.out"))
WT_100to250ns <- data.table::fread(text = readLines("3_WT/3_production/slurm-3548241.out")[sel])[, -1]
colnames(WT_100to250ns) <- columns
WT_100to250ns$TS <- WT_100to250ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_100to250ns <- WT_100to250ns[WT_100to250ns$TS %in% seq(WT_100to250ns$TS[1], WT_100to250ns$TS[length(WT_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("3_WT/3_production/slurm-7009945.out"))
WT_250to500ns <- data.table::fread(text = readLines("3_WT/3_production/slurm-7009945.out")[sel])[, -1]
colnames(WT_250to500ns) <- columns
WT_250to500ns$TS <- WT_250to500ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_250to500ns <- WT_250to500ns[WT_250to500ns$TS %in% seq(WT_250to500ns$TS[1], WT_250to500ns$TS[length(WT_250to500ns$TS)], 100000), ]

WT_eq_500ns <- rbind(WT_eq1, WT_eq2, WT_50ns, WT_50to100ns,
                     WT_100to250ns, WT_250to500ns)
WT_eq_500ns$time <- (WT_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(WT_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(WT_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("3_WT/4_production_repeat2/seokWT_100ns_r2.log"))
WT_r2_100ns <- data.table::fread(text = readLines("3_WT/4_production_repeat2/seokWT_100ns_r2.log")[sel])[, -1]
colnames(WT_r2_100ns) <- columns
WT_r2_100ns$TS <- WT_r2_100ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_r2_100ns <- WT_r2_100ns[WT_r2_100ns$TS %in% seq(WT_r2_100ns$TS[1], WT_r2_100ns$TS[length(WT_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("3_WT/4_production_repeat2/slurm-3548243.out"))
WT_r2_100to250ns <- data.table::fread(text = readLines("3_WT/4_production_repeat2/slurm-3548243.out")[sel])[, -1]
colnames(WT_r2_100to250ns) <- columns
WT_r2_100to250ns$TS <- WT_r2_100to250ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_r2_100to250ns <- WT_r2_100to250ns[WT_r2_100to250ns$TS %in% seq(WT_r2_100to250ns$TS[1], WT_r2_100to250ns$TS[length(WT_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("3_WT/4_production_repeat2/slurm-7009947.out"))
WT_r2_250to500ns <- data.table::fread(text = readLines("3_WT/4_production_repeat2/slurm-7009947.out")[sel])[, -1]
colnames(WT_r2_250to500ns) <- columns
WT_r2_250to500ns$TS <- WT_r2_250to500ns$TS + WT_eq2$TS[length(WT_eq2$TS)]
WT_r2_250to500ns <- WT_r2_250to500ns[WT_r2_250to500ns$TS %in% seq(WT_r2_250to500ns$TS[1], WT_r2_250to500ns$TS[length(WT_r2_250to500ns$TS)], 100000), ]

WT_r2_eq_500ns <- rbind(WT_eq1, WT_eq2, WT_r2_100ns,
                     WT_r2_100to250ns, WT_r2_250to500ns)
WT_r2_eq_500ns$time <- (WT_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(WT_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_wt <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                          top = text_grob("Wild-type", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "columns")))

# S50P --------------------------------------
sel <- grep("^ENERGY", readLines("7_S50P/1_eq1/seokS50P_eq1.log"))
S50P_eq1 <- data.table::fread(text = readLines("7_S50P/1_eq1/seokS50P_eq1.log")[sel])[, -1]
colnames(S50P_eq1) <- columns

sel <- grep("^ENERGY", readLines("7_S50P/2_eq2/seokS50P_eq2.log"))
S50P_eq2 <- data.table::fread(text = readLines("7_S50P/2_eq2/seokS50P_eq2.log")[sel])[, -1]
colnames(S50P_eq2) <- columns

sel <- grep("^ENERGY", readLines("7_S50P/3_production/seokS50P_100ns.log"))
S50P_100ns <- data.table::fread(text = readLines("7_S50P/3_production/seokS50P_100ns.log")[sel])[, -1]
colnames(S50P_100ns) <- columns
S50P_100ns$TS <- S50P_100ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_100ns <- S50P_100ns[S50P_100ns$TS %in% seq(S50P_100ns$TS[1], S50P_100ns$TS[length(S50P_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("7_S50P/3_production/slurm-3548255.out"))
S50P_100to250ns <- data.table::fread(text = readLines("7_S50P/3_production/slurm-3548255.out")[sel])[, -1]
colnames(S50P_100to250ns) <- columns
S50P_100to250ns$TS <- S50P_100to250ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_100to250ns <- S50P_100to250ns[S50P_100to250ns$TS %in% seq(S50P_100to250ns$TS[1], S50P_100to250ns$TS[length(S50P_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("7_S50P/3_production/slurm-7012220.out"))
S50P_250to500ns <- data.table::fread(text = readLines("7_S50P/3_production/slurm-7012220.out")[sel])[, -1]
colnames(S50P_250to500ns) <- columns
S50P_250to500ns$TS <- S50P_250to500ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_250to500ns <- S50P_250to500ns[S50P_250to500ns$TS %in% seq(S50P_250to500ns$TS[1], S50P_250to500ns$TS[length(S50P_250to500ns$TS)], 100000), ]

S50P_eq_500ns <- rbind(S50P_eq1, S50P_eq2, S50P_100ns,
                     S50P_100to250ns, S50P_250to500ns)
S50P_eq_500ns$time <- (S50P_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(S50P_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(S50P_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("7_S50P/4_production_repeat2/slurm-3290231.out"))
S50P_r2_100ns <- data.table::fread(text = readLines("7_S50P/4_production_repeat2/slurm-3290231.out")[sel])[, -1]
colnames(S50P_r2_100ns) <- columns
S50P_r2_100ns$TS <- S50P_r2_100ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_r2_100ns <- S50P_r2_100ns[S50P_r2_100ns$TS %in% seq(S50P_r2_100ns$TS[1], S50P_r2_100ns$TS[length(S50P_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("7_S50P/4_production_repeat2/slurm-3548257.out"))
S50P_r2_100to250ns <- data.table::fread(text = readLines("7_S50P/4_production_repeat2/slurm-3548257.out")[sel])[, -1]
colnames(S50P_r2_100to250ns) <- columns
S50P_r2_100to250ns$TS <- S50P_r2_100to250ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_r2_100to250ns <- S50P_r2_100to250ns[S50P_r2_100to250ns$TS %in% seq(S50P_r2_100to250ns$TS[1], S50P_r2_100to250ns$TS[length(S50P_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("7_S50P/4_production_repeat2/slurm-7012222.out"))
S50P_r2_250to500ns <- data.table::fread(text = readLines("7_S50P/4_production_repeat2/slurm-7012222.out")[sel])[, -1]
colnames(S50P_r2_250to500ns) <- columns
S50P_r2_250to500ns$TS <- S50P_r2_250to500ns$TS + S50P_eq2$TS[length(S50P_eq2$TS)]
S50P_r2_250to500ns <- S50P_r2_250to500ns[S50P_r2_250to500ns$TS %in% seq(S50P_r2_250to500ns$TS[1], S50P_r2_250to500ns$TS[length(S50P_r2_250to500ns$TS)], 100000), ]

S50P_r2_eq_500ns <- rbind(S50P_eq1, S50P_eq2, S50P_r2_100ns,
                        S50P_r2_100to250ns, S50P_r2_250to500ns)
S50P_r2_eq_500ns$time <- (S50P_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(S50P_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_S50P <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                          top = text_grob("p.Ser50Pro", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "fig_S50P", "columns")))

# R266Q --------------------------------------
sel <- grep("^ENERGY", readLines("8_R266Q/1_eq1/slurm-3277667.out"))
R266Q_eq1 <- data.table::fread(text = readLines("8_R266Q/1_eq1/slurm-3277667.out")[sel])[, -1]
colnames(R266Q_eq1) <- columns

sel <- grep("^ENERGY", readLines("8_R266Q/2_eq2/slurm-3282716.out"))
R266Q_eq2 <- data.table::fread(text = readLines("8_R266Q/2_eq2/slurm-3282716.out")[sel])[, -1]
colnames(R266Q_eq2) <- columns

sel <- grep("^ENERGY", readLines("8_R266Q/3_production/slurm-3287397.out"))
R266Q_100ns <- data.table::fread(text = readLines("8_R266Q/3_production/slurm-3287397.out")[sel])[, -1]
colnames(R266Q_100ns) <- columns
R266Q_100ns$TS <- R266Q_100ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_100ns <- R266Q_100ns[R266Q_100ns$TS %in% seq(R266Q_100ns$TS[1], R266Q_100ns$TS[length(R266Q_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("8_R266Q/3_production/slurm-3548260.out"))
R266Q_100to250ns <- data.table::fread(text = readLines("8_R266Q/3_production/slurm-3548260.out")[sel])[, -1]
colnames(R266Q_100to250ns) <- columns
R266Q_100to250ns$TS <- R266Q_100to250ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_100to250ns <- R266Q_100to250ns[R266Q_100to250ns$TS %in% seq(R266Q_100to250ns$TS[1], R266Q_100to250ns$TS[length(R266Q_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("8_R266Q/3_production/slurm-7012223.out"))
R266Q_250to500ns <- data.table::fread(text = readLines("8_R266Q/3_production/slurm-7012223.out")[sel])[, -1]
colnames(R266Q_250to500ns) <- columns
R266Q_250to500ns$TS <- R266Q_250to500ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_250to500ns <- R266Q_250to500ns[R266Q_250to500ns$TS %in% seq(R266Q_250to500ns$TS[1], R266Q_250to500ns$TS[length(R266Q_250to500ns$TS)], 100000), ]

R266Q_eq_500ns <- rbind(R266Q_eq1, R266Q_eq2, R266Q_100ns,
                       R266Q_100to250ns, R266Q_250to500ns)
R266Q_eq_500ns$time <- (R266Q_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(R266Q_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(R266Q_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("8_R266Q/4_production_repeat2/slurm-3287398.out"))
R266Q_r2_100ns <- data.table::fread(text = readLines("8_R266Q/4_production_repeat2/slurm-3287398.out")[sel])[, -1]
colnames(R266Q_r2_100ns) <- columns
R266Q_r2_100ns$TS <- R266Q_r2_100ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_r2_100ns <- R266Q_r2_100ns[R266Q_r2_100ns$TS %in% seq(R266Q_r2_100ns$TS[1], R266Q_r2_100ns$TS[length(R266Q_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("8_R266Q/4_production_repeat2/slurm-3548261.out"))
R266Q_r2_100to250ns <- data.table::fread(text = readLines("8_R266Q/4_production_repeat2/slurm-3548261.out")[sel])[, -1]
colnames(R266Q_r2_100to250ns) <- columns
R266Q_r2_100to250ns$TS <- R266Q_r2_100to250ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_r2_100to250ns <- R266Q_r2_100to250ns[R266Q_r2_100to250ns$TS %in% seq(R266Q_r2_100to250ns$TS[1], R266Q_r2_100to250ns$TS[length(R266Q_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("8_R266Q/4_production_repeat2/slurm-7012225.out"))
R266Q_r2_250to500ns <- data.table::fread(text = readLines("8_R266Q/4_production_repeat2/slurm-7012225.out")[sel])[, -1]
colnames(R266Q_r2_250to500ns) <- columns
R266Q_r2_250to500ns$TS <- R266Q_r2_250to500ns$TS + R266Q_eq2$TS[length(R266Q_eq2$TS)]
R266Q_r2_250to500ns <- R266Q_r2_250to500ns[R266Q_r2_250to500ns$TS %in% seq(R266Q_r2_250to500ns$TS[1], R266Q_r2_250to500ns$TS[length(R266Q_r2_250to500ns$TS)], 100000), ]

R266Q_r2_eq_500ns <- rbind(R266Q_eq1, R266Q_eq2, R266Q_r2_100ns,
                          R266Q_r2_100to250ns, R266Q_r2_250to500ns)
R266Q_r2_eq_500ns$time <- (R266Q_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(R266Q_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_R266Q <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                            top = text_grob("p.Arg266Gln", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "fig_S50P", "fig_R266Q", "columns")))

# G317R --------------------------------------
sel <- grep("^ENERGY", readLines("5_G317R/1_eq1/seokG317R_eq1.log"))
G317R_eq1 <- data.table::fread(text = readLines("5_G317R/1_eq1/seokG317R_eq1.log")[sel])[, -1]
colnames(G317R_eq1) <- columns

sel <- grep("^ENERGY", readLines("5_G317R/2_eq2/seokG317R_eq2.log"))
G317R_eq2 <- data.table::fread(text = readLines("5_G317R/2_eq2/seokG317R_eq2.log")[sel])[, -1]
colnames(G317R_eq2) <- columns

sel <- grep("^ENERGY", readLines("5_G317R/3_production/seokG317R_100ns.log"))
G317R_100ns <- data.table::fread(text = readLines("5_G317R/3_production/seokG317R_100ns.log")[sel])[, -1]
colnames(G317R_100ns) <- columns
G317R_100ns$TS <- G317R_100ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_100ns <- G317R_100ns[G317R_100ns$TS %in% seq(G317R_100ns$TS[1], G317R_100ns$TS[length(G317R_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("5_G317R/3_production/slurm-3551995.out"))
G317R_100to250ns <- data.table::fread(text = readLines("5_G317R/3_production/slurm-3551995.out")[sel])[, -1]
colnames(G317R_100to250ns) <- columns
G317R_100to250ns$TS <- G317R_100to250ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_100to250ns <- G317R_100to250ns[G317R_100to250ns$TS %in% seq(G317R_100to250ns$TS[1], G317R_100to250ns$TS[length(G317R_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("5_G317R/3_production/slurm-7010497.out"))
G317R_250to500ns <- data.table::fread(text = readLines("5_G317R/3_production/slurm-7010497.out")[sel])[, -1]
colnames(G317R_250to500ns) <- columns
G317R_250to500ns$TS <- G317R_250to500ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_250to500ns <- G317R_250to500ns[G317R_250to500ns$TS %in% seq(G317R_250to500ns$TS[1], G317R_250to500ns$TS[length(G317R_250to500ns$TS)], 100000), ]

G317R_eq_500ns <- rbind(G317R_eq1, G317R_eq2, G317R_100ns,
                        G317R_100to250ns, G317R_250to500ns)
G317R_eq_500ns$time <- (G317R_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(G317R_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(G317R_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("5_G317R/4_production_repeat2/slurm-3290283.out"))
G317R_r2_100ns <- data.table::fread(text = readLines("5_G317R/4_production_repeat2/slurm-3290283.out")[sel])[, -1]
colnames(G317R_r2_100ns) <- columns
G317R_r2_100ns$TS <- G317R_r2_100ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_r2_100ns <- G317R_r2_100ns[G317R_r2_100ns$TS %in% seq(G317R_r2_100ns$TS[1], G317R_r2_100ns$TS[length(G317R_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("5_G317R/4_production_repeat2/slurm-3551996.out"))
G317R_r2_100to250ns <- data.table::fread(text = readLines("5_G317R/4_production_repeat2/slurm-3551996.out")[sel])[, -1]
colnames(G317R_r2_100to250ns) <- columns
G317R_r2_100to250ns$TS <- G317R_r2_100to250ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_r2_100to250ns <- G317R_r2_100to250ns[G317R_r2_100to250ns$TS %in% seq(G317R_r2_100to250ns$TS[1], G317R_r2_100to250ns$TS[length(G317R_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("5_G317R/4_production_repeat2/slurm-7010498.out"))
G317R_r2_250to500ns <- data.table::fread(text = readLines("5_G317R/4_production_repeat2/slurm-7010498.out")[sel])[, -1]
colnames(G317R_r2_250to500ns) <- columns
G317R_r2_250to500ns$TS <- G317R_r2_250to500ns$TS + G317R_eq2$TS[length(G317R_eq2$TS)]
G317R_r2_250to500ns <- G317R_r2_250to500ns[G317R_r2_250to500ns$TS %in% seq(G317R_r2_250to500ns$TS[1], G317R_r2_250to500ns$TS[length(G317R_r2_250to500ns$TS)], 100000), ]

G317R_r2_eq_500ns <- rbind(G317R_eq1, G317R_eq2, G317R_r2_100ns,
                           G317R_r2_100to250ns, G317R_r2_250to500ns)
G317R_r2_eq_500ns$time <- (G317R_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(G317R_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_G317R <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                             top = text_grob("p.Gly317Arg", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "fig_S50P", "fig_R266Q",
                          "fig_G317R", "columns")))

# P389L --------------------------------------
sel <- grep("^ENERGY", readLines("4_P389L/1_eq1/seokP389L_eq1.log"))
P389L_eq1 <- data.table::fread(text = readLines("4_P389L/1_eq1/seokP389L_eq1.log")[sel])[, -1]
colnames(P389L_eq1) <- columns

sel <- grep("^ENERGY", readLines("4_P389L/2_eq2/seokP389L_eq2.log"))
P389L_eq2 <- data.table::fread(text = readLines("4_P389L/2_eq2/seokP389L_eq2.log")[sel])[, -1]
colnames(P389L_eq2) <- columns

sel <- grep("^ENERGY", readLines("4_P389L/3_production/seokP389L_50ns.log"))
P389L_50ns <- data.table::fread(text = readLines("4_P389L/3_production/seokP389L_50ns.log")[sel])[, -1]
colnames(P389L_50ns) <- columns
P389L_50ns$TS <- P389L_50ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_50ns <- P389L_50ns[P389L_50ns$TS %in% seq(P389L_50ns$TS[1], P389L_50ns$TS[length(P389L_50ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("4_P389L/3_production/seokP389L_50to100ns.log"))
P389L_50to100ns <- data.table::fread(text = readLines("4_P389L/3_production/seokP389L_50to100ns.log")[sel])[, -1]
colnames(P389L_50to100ns) <- columns
P389L_50to100ns$TS <- P389L_50to100ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_50to100ns <- P389L_50to100ns[P389L_50to100ns$TS %in% seq(P389L_50to100ns$TS[1], P389L_50to100ns$TS[length(P389L_50to100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("4_P389L/3_production/slurm-3548245.out"))
P389L_100to250ns <- data.table::fread(text = readLines("4_P389L/3_production/slurm-3548245.out")[sel])[, -1]
colnames(P389L_100to250ns) <- columns
P389L_100to250ns$TS <- P389L_100to250ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_100to250ns <- P389L_100to250ns[P389L_100to250ns$TS %in% seq(P389L_100to250ns$TS[1], P389L_100to250ns$TS[length(P389L_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("4_P389L/3_production/slurm-7009955.out"))
P389L_250to500ns <- data.table::fread(text = readLines("4_P389L/3_production/slurm-7009955.out")[sel])[, -1]
colnames(P389L_250to500ns) <- columns
P389L_250to500ns$TS <- P389L_250to500ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_250to500ns <- P389L_250to500ns[P389L_250to500ns$TS %in% seq(P389L_250to500ns$TS[1], P389L_250to500ns$TS[length(P389L_250to500ns$TS)], 100000), ]

P389L_eq_500ns <- rbind(P389L_eq1, P389L_eq2, P389L_50ns, P389L_50to100ns,
                        P389L_100to250ns, P389L_250to500ns)
P389L_eq_500ns$time <- (P389L_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(P389L_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(P389L_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("4_P389L/4_production_repeat2/slurm-3290296.out"))
P389L_r2_100ns <- data.table::fread(text = readLines("4_P389L/4_production_repeat2/slurm-3290296.out")[sel])[, -1]
colnames(P389L_r2_100ns) <- columns
P389L_r2_100ns$TS <- P389L_r2_100ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_r2_100ns <- P389L_r2_100ns[P389L_r2_100ns$TS %in% seq(P389L_r2_100ns$TS[1], P389L_r2_100ns$TS[length(P389L_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("4_P389L/4_production_repeat2/slurm-3548246.out"))
P389L_r2_100to250ns <- data.table::fread(text = readLines("4_P389L/4_production_repeat2/slurm-3548246.out")[sel])[, -1]
colnames(P389L_r2_100to250ns) <- columns
P389L_r2_100to250ns$TS <- P389L_r2_100to250ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_r2_100to250ns <- P389L_r2_100to250ns[P389L_r2_100to250ns$TS %in% seq(P389L_r2_100to250ns$TS[1], P389L_r2_100to250ns$TS[length(P389L_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("4_P389L/4_production_repeat2/slurm-7009959.out"))
P389L_r2_250to500ns <- data.table::fread(text = readLines("4_P389L/4_production_repeat2/slurm-7009959.out")[sel])[, -1]
colnames(P389L_r2_250to500ns) <- columns
P389L_r2_250to500ns$TS <- P389L_r2_250to500ns$TS + P389L_eq2$TS[length(P389L_eq2$TS)]
P389L_r2_250to500ns <- P389L_r2_250to500ns[P389L_r2_250to500ns$TS %in% seq(P389L_r2_250to500ns$TS[1], P389L_r2_250to500ns$TS[length(P389L_r2_250to500ns$TS)], 100000), ]

P389L_r2_eq_500ns <- rbind(P389L_eq1, P389L_eq2, P389L_r2_100ns,
                           P389L_r2_100to250ns, P389L_r2_250to500ns)
P389L_r2_eq_500ns$time <- (P389L_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(P389L_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_P389L <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                             top = text_grob("p.Pro389Leu", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "fig_S50P", "fig_R266Q",
                          "fig_G317R", "fig_P389L", "columns")))

# L430P --------------------------------------
sel <- grep("^ENERGY", readLines("6_L430P/1_eq1/seokL430P_eq1.log"))
L430P_eq1 <- data.table::fread(text = readLines("6_L430P/1_eq1/seokL430P_eq1.log")[sel])[, -1]
colnames(L430P_eq1) <- columns

sel <- grep("^ENERGY", readLines("6_L430P/2_eq2/seokL430P_eq2.log"))
L430P_eq2 <- data.table::fread(text = readLines("6_L430P/2_eq2/seokL430P_eq2.log")[sel])[, -1]
colnames(L430P_eq2) <- columns

sel <- grep("^ENERGY", readLines("6_L430P/3_production/seokL430P_100ns.log"))
L430P_100ns <- data.table::fread(text = readLines("6_L430P/3_production/seokL430P_100ns.log")[sel])[, -1]
colnames(L430P_100ns) <- columns
L430P_100ns$TS <- L430P_100ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_100ns <- L430P_100ns[L430P_100ns$TS %in% seq(L430P_100ns$TS[1], L430P_100ns$TS[length(L430P_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("6_L430P/3_production/slurm-3551998.out"))
L430P_100to250ns <- data.table::fread(text = readLines("6_L430P/3_production/slurm-3551998.out")[sel])[, -1]
colnames(L430P_100to250ns) <- columns
L430P_100to250ns$TS <- L430P_100to250ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_100to250ns <- L430P_100to250ns[L430P_100to250ns$TS %in% seq(L430P_100to250ns$TS[1], L430P_100to250ns$TS[length(L430P_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("6_L430P/3_production/slurm-7012218.out"))
L430P_250to500ns <- data.table::fread(text = readLines("6_L430P/3_production/slurm-7012218.out")[sel])[, -1]
colnames(L430P_250to500ns) <- columns
L430P_250to500ns$TS <- L430P_250to500ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_250to500ns <- L430P_250to500ns[L430P_250to500ns$TS %in% seq(L430P_250to500ns$TS[1], L430P_250to500ns$TS[length(L430P_250to500ns$TS)], 100000), ]

L430P_eq_500ns <- rbind(L430P_eq1, L430P_eq2, L430P_100ns,
                        L430P_100to250ns, L430P_250to500ns)
L430P_eq_500ns$time <- (L430P_eq_500ns$TS * 2) / 1e+06

p1 <- ggplot(L430P_eq_500ns, aes(x = time, y = POTENTIAL))
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p1 <- p1 + ggtitle("Minimization and Equilibration")
p1 <- p1 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1 <- gg.gap(plot = p1,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

p2 <- ggplot(L430P_eq_500ns, aes(x = time, y = POTENTIAL))
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p2 <- p2 + ggtitle("Repeat 1")
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2 <- gg.gap(plot = p2,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))


sel <- grep("^ENERGY", readLines("6_L430P/4_production_repeat2/slurm-3290254.out"))
L430P_r2_100ns <- data.table::fread(text = readLines("6_L430P/4_production_repeat2/slurm-3290254.out")[sel])[, -1]
colnames(L430P_r2_100ns) <- columns
L430P_r2_100ns$TS <- L430P_r2_100ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_r2_100ns <- L430P_r2_100ns[L430P_r2_100ns$TS %in% seq(L430P_r2_100ns$TS[1], L430P_r2_100ns$TS[length(L430P_r2_100ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("6_L430P/4_production_repeat2/slurm-3552000.out"))
L430P_r2_100to250ns <- data.table::fread(text = readLines("6_L430P/4_production_repeat2/slurm-3552000.out")[sel])[, -1]
colnames(L430P_r2_100to250ns) <- columns
L430P_r2_100to250ns$TS <- L430P_r2_100to250ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_r2_100to250ns <- L430P_r2_100to250ns[L430P_r2_100to250ns$TS %in% seq(L430P_r2_100to250ns$TS[1], L430P_r2_100to250ns$TS[length(L430P_r2_100to250ns$TS)], 100000), ]

sel <- grep("^ENERGY", readLines("6_L430P/4_production_repeat2/slurm-7012219.out"))
L430P_r2_250to500ns <- data.table::fread(text = readLines("6_L430P/4_production_repeat2/slurm-7012219.out")[sel])[, -1]
colnames(L430P_r2_250to500ns) <- columns
L430P_r2_250to500ns$TS <- L430P_r2_250to500ns$TS + L430P_eq2$TS[length(L430P_eq2$TS)]
L430P_r2_250to500ns <- L430P_r2_250to500ns[L430P_r2_250to500ns$TS %in% seq(L430P_r2_250to500ns$TS[1], L430P_r2_250to500ns$TS[length(L430P_r2_250to500ns$TS)], 100000), ]

L430P_r2_eq_500ns <- rbind(L430P_eq1, L430P_eq2, L430P_r2_100ns,
                           L430P_r2_100to250ns, L430P_r2_250to500ns)
L430P_r2_eq_500ns$time <- (L430P_r2_eq_500ns$TS * 2) / 1e+06

p3 <- ggplot(L430P_r2_eq_500ns, aes(x = time, y = POTENTIAL))
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + labs(x = "Time (ns)", y = "Potential Energy (kcal/mol)")
p3 <- p3 + ggtitle("Repeat 2")
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 10, 0.25), limits = c(0, 2))
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_text(size = 12),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3 <- gg.gap(plot = p3,
             segments = c(-87000, 750000),
             ylim = c(-625000, 50925000),
             tick_width = c(100000, 1e+07))

fig_L430P <- annotate_figure(ggarrange(p1, p2, p3, ncol = 1, nrow = 3),
                             top = text_grob("p.Leu430Pro", face = "bold", size = 24))

rm(list = setdiff(ls(), c("fig_wt", "fig_S50P", "fig_R266Q",
                          "fig_G317R", "fig_P389L", "fig_L430P", "columns")))



# out ------------------
pdf("_outputs/Potential_Energy_500ns.pdf", width = 12, height = 14)
fig_wt
fig_S50P
fig_R266Q
fig_G317R
fig_P389L
fig_L430P
dev.off()
