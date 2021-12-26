##################################################
## Project: M1AP mutations
## Script purpose: visualization of analyses
## Date: Dec 22, 2021
## Author: Umut Gerlevik
##################################################

# Read libraries -----------
library(ggplot2)
library(ggpubr)

# Set working directory -----------
setwd("D:/SezermanLab/M1AP/2_simulations/GalaxyWeb_Seok/")

# Common variables ----------
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")
colorPalet <- grDevices::colors()[grDevices::colors() %in% c("steelblue3", "royalblue2",
                                                             "orangered3", "indianred3")]

# RMSD ------------
rmsd_wt_r1 <- read.table("3_WT/3_production/rmsd_seokWT.dat", sep = ",")
rmsd_wt_r2 <- read.table("3_WT/4_production_repeat2/rmsd_seokWT.dat", sep = ",")
rmsd_s50p_r1 <- read.table("7_S50P/3_production/rmsd_seokS50P.dat", sep = ",")
rmsd_s50p_r2 <- read.table("7_S50P/4_production_repeat2/rmsd_seokS50P.dat", sep = ",")
rmsd_r266q_r1 <- read.table("8_R266Q/3_production/rmsd_seokR266Q.dat", sep = ",")
rmsd_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/rmsd_seokR266Q.dat", sep = ",")
rmsd_g317r_r1 <- read.table("5_G317R/3_production/rmsd_seokG317R.dat", sep = ",")
rmsd_g317r_r2 <- read.table("5_G317R/4_production_repeat2/rmsd_seokG317R.dat", sep = ",")
rmsd_p389l_r1 <- read.table("4_P389L/3_production/rmsd_seokP389L.dat", sep = ",")
rmsd_p389l_r2 <- read.table("4_P389L/4_production_repeat2/rmsd_seokP389L.dat", sep = ",")
rmsd_l430p_r1 <- read.table("6_L430P/3_production/rmsd_seokL430P.dat", sep = ",")
rmsd_l430p_r2 <- read.table("6_L430P/4_production_repeat2/rmsd_seokL430P.dat", sep = ",")

rmsd_r1 <- cbind(rmsd_wt_r1, rmsd_s50p_r1[2], rmsd_r266q_r1[2], rmsd_g317r_r1[2],
                 rmsd_p389l_r1[2], rmsd_l430p_r1[2])
rmsd_r1[1] <- seq(0, 499.99, 0.1)
colnames(rmsd_r1) <- c("time", "Wild-type_R1", "p.Ser50Pro_R1", "p.Arg266Gln_R1", "p.Gly317Arg_R1",
                       "p.Pro389Leu_R1", "p.Leu430Pro_R1")
rmsd_r1 <- reshape2::melt(rmsd_r1[1:7], id.var = "time")
rmsd_r1$rep <- "Repeat 1"

rmsd_r2 <- cbind(rmsd_wt_r2, rmsd_s50p_r2[2], rmsd_r266q_r2[2], rmsd_g317r_r2[2],
                 rmsd_p389l_r2[2], rmsd_l430p_r2[2])
rmsd_r2[1] <- seq(0, 499.99, 0.1)
colnames(rmsd_r2) <- c("time", "Wild-type_R2", "p.Ser50Pro_R2", "p.Arg266Gln_R2", "p.Gly317Arg_R2",
                       "p.Pro389Leu_R2", "p.Leu430Pro_R2")
rmsd_r2 <- reshape2::melt(rmsd_r2[1:7], id.var = "time")
rmsd_r2$rep <- "Repeat 2"

rmsd_all <- rbind(rmsd_r1, rmsd_r2)
rmsd_all$variable <- factor(rmsd_all$variable, levels = c("Wild-type_R1", "Wild-type_R2", "p.Ser50Pro_R1","p.Ser50Pro_R2",
                                                          "p.Arg266Gln_R1", "p.Arg266Gln_R2", "p.Gly317Arg_R1", "p.Gly317Arg_R2",
                                                          "p.Pro389Leu_R1", "p.Pro389Leu_R2", "p.Leu430Pro_R1", "p.Leu430Pro_R2"))

p1 <- ggplot(rmsd_all[rmsd_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                             paste0(mutationlist[1], "_R1"),
                                             paste0(mutationlist[1], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 500, 50))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 50, 1), limits = c(0, 15))
p1 <- p1 + labs(x = "Time (ns)", y = "RMSD (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet, breaks = c(paste0(mutationlist[1], "_R1"),
                                                              paste0(mutationlist[1], "_R2"),
                                                              "Wild-type_R1", "Wild-type_R2"))
# p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6.5),
                 legend.key.height = unit(0.3, 'cm'),
                 legend.position = c(0.89, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(rmsd_all[rmsd_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                             paste0(mutationlist[2], "_R1"),
                                             paste0(mutationlist[2], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 500, 50))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 50, 1), limits = c(0, 15))
p2 <- p2 + labs(x = "Time (ns)", y = "RMSD (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet, breaks = c(paste0(mutationlist[2], "_R1"),
                                                              paste0(mutationlist[2], "_R2"),
                                                              "Wild-type_R1", "Wild-type_R2"))
# p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6.5),
                 legend.key.height = unit(0.3, 'cm'),
                 legend.position = c(0.89, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(rmsd_all[rmsd_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                             paste0(mutationlist[3], "_R1"),
                                             paste0(mutationlist[3], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 500, 50))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 50, 1), limits = c(0, 15))
p3 <- p3 + labs(x = "Time (ns)", y = "RMSD (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet, breaks = c(paste0(mutationlist[3], "_R1"),
                                                              paste0(mutationlist[3], "_R2"),
                                                              "Wild-type_R1", "Wild-type_R2"))
# p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6.5),
                 legend.key.height = unit(0.3, 'cm'),
                 legend.position = c(0.89, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(rmsd_all[rmsd_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[4], "_R1"),
                                           paste0(mutationlist[4], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7, alpha = 0.75) 
p4 <- p4 + geom_line(size = 0.65, alpha = 0.75)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 500, 50))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 50, 1), limits = c(0, 15))
p4 <- p4 + labs(x = "Time (ns)", y = "RMSD (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet, breaks = c(paste0(mutationlist[4], "_R1"),
                                                              paste0(mutationlist[4], "_R2"),
                                                              "Wild-type_R1", "Wild-type_R2"))
# p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6.5),
                 legend.key.height = unit(0.3, 'cm'),
                 legend.position = c(0.89, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(rmsd_all[rmsd_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                             paste0(mutationlist[5], "_R1"),
                                             paste0(mutationlist[5], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7, alpha = 0.75) 
p5 <- p5 + geom_line(size = 0.65, alpha = 0.75)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 500, 50))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 50, 1), limits = c(0, 15))
p5 <- p5 + labs(x = "Time (ns)", y = "RMSD (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet, breaks = c(paste0(mutationlist[5], "_R1"),
                                                              paste0(mutationlist[5], "_R2"),
                                                              "Wild-type_R1", "Wild-type_R2"))
# p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 9),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6.5),
                 legend.key.height = unit(0.3, 'cm'),
                 legend.position = c(0.89, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/RMSD_500ns.pdf", width = 6, height = 12)
annotate_figure(fig, left = text_grob("Root-mean-square Deviation (Å)", face = "bold", size = 12, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 12))
dev.off()

