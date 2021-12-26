##################################################
## Project: M1AP mutations
## Script purpose: Visualization of H-bonds
## Date: Dec 24, 2021
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

# hbonds --------------
hbonds_wt_r1 <- read.table("3_WT/3_production/hbonds_seokWT.dat")
hbonds_wt_r2 <- read.table("3_WT/4_production_repeat2/hbonds_seokWT.dat")
hbonds_s50p_r1 <- read.table("7_S50P/3_production/hbonds_seokS50P.dat")
hbonds_s50p_r2 <- read.table("7_S50P/4_production_repeat2/hbonds_seokS50P.dat")
hbonds_r266q_r1 <- read.table("8_R266Q/3_production/hbonds_seokR266Q.dat")
hbonds_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/hbonds_seokR266Q.dat")
hbonds_g317r_r1 <- read.table("5_G317R/3_production/hbonds_seokG317R.dat")
hbonds_g317r_r2 <- read.table("5_G317R/4_production_repeat2/hbonds_seokG317R.dat")
hbonds_p389l_r1 <- read.table("4_P389L/3_production/hbonds_seokP389L.dat")
hbonds_p389l_r2 <- read.table("4_P389L/4_production_repeat2/hbonds_seokP389L.dat")
hbonds_l430p_r1 <- read.table("6_L430P/3_production/hbonds_seokL430P.dat")
hbonds_l430p_r2 <- read.table("6_L430P/4_production_repeat2/hbonds_seokL430P.dat")

hbonds_r1 <- cbind(hbonds_wt_r1, hbonds_s50p_r1[2], hbonds_r266q_r1[2], hbonds_g317r_r1[2],
                hbonds_p389l_r1[2], hbonds_l430p_r1[2])
hbonds_r1[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_r1) <- c("time", "Wild-type_R1", "p.Ser50Pro_R1", "p.Arg266Gln_R1", "p.Gly317Arg_R1",
                      "p.Pro389Leu_R1", "p.Leu430Pro_R1")
hbonds_r1 <- reshape2::melt(hbonds_r1[1:7], id.var = "time")
hbonds_r1$rep <- "Repeat 1"

hbonds_r2 <- cbind(hbonds_wt_r2, hbonds_s50p_r2[2], hbonds_r266q_r2[2], hbonds_g317r_r2[2],
                hbonds_p389l_r2[2], hbonds_l430p_r2[2])
hbonds_r2[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_r2) <- c("time", "Wild-type_R2", "p.Ser50Pro_R2", "p.Arg266Gln_R2", "p.Gly317Arg_R2",
                      "p.Pro389Leu_R2", "p.Leu430Pro_R2")
hbonds_r2 <- reshape2::melt(hbonds_r2[1:7], id.var = "time")
hbonds_r2 $rep <- "Repeat 2"

hbonds_all <- rbind(hbonds_r1, hbonds_r2)

hbonds_all$variable <- factor(hbonds_all$variable, levels = c("Wild-type_R1", "Wild-type_R2", "p.Ser50Pro_R1","p.Ser50Pro_R2",
                                                        "p.Arg266Gln_R1", "p.Arg266Gln_R2", "p.Gly317Arg_R1", "p.Gly317Arg_R2",
                                                        "p.Pro389Leu_R1", "p.Pro389Leu_R2", "p.Leu430Pro_R1", "p.Leu430Pro_R2"))

p1 <- ggplot(hbonds_all[hbonds_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[1], "_R1"),
                                           paste0(mutationlist[1], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7, alpha = 0.75) 
p1 <- p1 + geom_line(size = 0.65, alpha = 0.75)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 500, 50))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 500, 10), limits = c(60, 140))
p1 <- p1 + labs(x = "Time (ns)", y = "hbonds (Å)")
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
                 legend.position = c(0.11, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(hbonds_all[hbonds_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[2], "_R1"),
                                           paste0(mutationlist[2], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7, alpha = 0.75) 
p2 <- p2 + geom_line(size = 0.65, alpha = 0.75)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 500, 50))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 500, 10), limits = c(60, 140))
p2 <- p2 + labs(x = "Time (ns)", y = "hbonds (Å)")
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
                 legend.position = c(0.11, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(hbonds_all[hbonds_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[3], "_R1"),
                                           paste0(mutationlist[3], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7, alpha = 0.75) 
p3 <- p3 + geom_line(size = 0.65, alpha = 0.75)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 500, 50))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 500, 10), limits = c(60, 140))
p3 <- p3 + labs(x = "Time (ns)", y = "hbonds (Å)")
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
                 legend.position = c(0.11, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(hbonds_all[hbonds_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[4], "_R1"),
                                           paste0(mutationlist[4], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7, alpha = 0.75) 
p4 <- p4 + geom_line(size = 0.65, alpha = 0.75)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 500, 50))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 500, 10), limits = c(60, 140))
p4 <- p4 + labs(x = "Time (ns)", y = "hbonds (Å)")
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
                 legend.position = c(0.11, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(hbonds_all[hbonds_all$variable %in% c("Wild-type_R1", "Wild-type_R2",
                                           paste0(mutationlist[5], "_R1"),
                                           paste0(mutationlist[5], "_R2")), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7, alpha = 0.75) 
p5 <- p5 + geom_line(size = 0.65, alpha = 0.75)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 500, 50))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 500, 10), limits = c(60, 140))
p5 <- p5 + labs(x = "Time (ns)", y = "hbonds (Å)")
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
                 legend.position = c(0.11, 0.89),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/Hbonds_500ns.pdf", width = 6, height = 12)
annotate_figure(fig, left = text_grob("Number of H-bonds", face = "bold", size = 12, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 12))
dev.off()

