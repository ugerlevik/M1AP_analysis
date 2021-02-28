##################################################
## Project: M1AP mutations
## Script purpose: visualization of analyses
## Date: August 11, 2020
## Author: Umut Gerlevik
##################################################

# Read libraries -----------
library(ggplot2)
library(ggpubr)

# Common variables ----------
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")
palet <- ggsci::pal_jama()
colorPalet <- palet(6)

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
rmsd_r1[1] <- seq(0, 249.99, 0.1)
colnames(rmsd_r1) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
rmsd_r1 <- reshape2::melt(rmsd_r1[1:7], id.var = "time")
rmsd_r1$rep <- "Repeat 1"

rmsd_r2 <- cbind(rmsd_wt_r2, rmsd_s50p_r2[2], rmsd_r266q_r2[2], rmsd_g317r_r2[2],
                 rmsd_p389l_r2[2], rmsd_l430p_r2[2])
rmsd_r2[1] <- seq(0, 249.99, 0.1)
colnames(rmsd_r2) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
rmsd_r2 <- reshape2::melt(rmsd_r2[1:7], id.var = "time")
rmsd_r2$rep <- "Repeat 2"

rmsd_all <- rbind(rmsd_r1, rmsd_r2)

p1 <- ggplot(rmsd_all[rmsd_all$variable == c("Wild-type", mutationlist[1]), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 300, 25))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 50, 1))
p1 <- p1 + labs(x = "Time (ns)", y = "RMSD (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.9, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(rmsd_all[rmsd_all$variable == c("Wild-type", mutationlist[2]), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 300, 25))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 50, 1))
p2 <- p2 + labs(x = "Time (ns)", y = "RMSD (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.9, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(rmsd_all[rmsd_all$variable == c("Wild-type", mutationlist[3]), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 300, 25))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 50, 1))
p3 <- p3 + labs(x = "Time (ns)", y = "RMSD (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.9, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(rmsd_all[rmsd_all$variable == c("Wild-type", mutationlist[4]), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 300, 25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 50, 1))
p4 <- p4 + labs(x = "Time (ns)", y = "RMSD (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.9, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(rmsd_all[rmsd_all$variable == c("Wild-type", mutationlist[5]), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 300, 25))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 50, 1))
p5 <- p5 + labs(x = "Time (ns)", y = "RMSD (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.9, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/RMSD.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("Root-mean-square Deviation (Å)", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 20))
dev.off()

# RMSF -------------
rmsf_wt_r1 <- read.table("3_WT/3_production/rmsf_seokWT.dat", sep = ",")
rmsf_wt_r2 <- read.table("3_WT/4_production_repeat2/rmsf_seokWT.dat", sep = ",")
rmsf_s50p_r1 <- read.table("7_S50P/3_production/rmsf_seokS50P.dat", sep = ",")
rmsf_s50p_r2 <- read.table("7_S50P/4_production_repeat2/rmsf_seokS50P.dat", sep = ",")
rmsf_r266q_r1 <- read.table("8_R266Q/3_production/rmsf_seokR266Q.dat", sep = ",")
rmsf_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/rmsf_seokR266Q.dat", sep = ",")
rmsf_g317r_r1 <- read.table("5_G317R/3_production/rmsf_seokG317R.dat", sep = ",")
rmsf_g317r_r2 <- read.table("5_G317R/4_production_repeat2/rmsf_seokG317R.dat", sep = ",")
rmsf_p389l_r1 <- read.table("4_P389L/3_production/rmsf_seokP389L.dat", sep = ",")
rmsf_p389l_r2 <- read.table("4_P389L/4_production_repeat2/rmsf_seokP389L.dat", sep = ",")
rmsf_l430p_r1 <- read.table("6_L430P/3_production/rmsf_seokL430P.dat", sep = ",")
rmsf_l430p_r2 <- read.table("6_L430P/4_production_repeat2/rmsf_seokL430P.dat", sep = ",")

rmsf_r1 <- cbind(rmsf_wt_r1, rmsf_s50p_r1[2], rmsf_r266q_r1[2], rmsf_g317r_r1[2],
                 rmsf_p389l_r1[2], rmsf_l430p_r1[2])
colnames(rmsf_r1) <- c("resid", "Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
rmsf_r1 <- reshape2::melt(rmsf_r1[1:7], id.var = "resid")
rmsf_r1$rep <- "Repeat 1"

rmsf_r2 <- cbind(rmsf_wt_r2, rmsf_s50p_r2[2], rmsf_r266q_r2[2], rmsf_g317r_r2[2],
                 rmsf_p389l_r2[2], rmsf_l430p_r2[2])
colnames(rmsf_r2) <- c("resid", "Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
rmsf_r2 <- reshape2::melt(rmsf_r2[1:7], id.var = "resid")
rmsf_r2 $rep <- "Repeat 2"

rmsf_all <- rbind(rmsf_r1, rmsf_r2)

p1 <- ggplot(rmsf_all[rmsf_all$variable == c("Wild-type", mutationlist[1]), ],
             aes(x=resid, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = c(1, seq(50, 500, 50), 530))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 50, 2))
p1 <- p1 + labs(x = "Residue Index", y = "RMSF (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle = 30),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(rmsf_all[rmsf_all$variable == c("Wild-type", mutationlist[2]), ],
             aes(x=resid, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = c(1, seq(50, 500, 50), 530))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 50, 2))
p2 <- p2 + labs(x = "Residue Index", y = "RMSF (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle = 30),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(rmsf_all[rmsf_all$variable == c("Wild-type", mutationlist[3]), ],
             aes(x=resid, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = c(1, seq(50, 500, 50), 530))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 50, 2))
p3 <- p3 + labs(x = "Residue Index", y = "RMSF (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle = 30),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(rmsf_all[rmsf_all$variable == c("Wild-type", mutationlist[4]), ],
             aes(x=resid, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = c(1, seq(50, 500, 50), 530))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 50, 2))
p4 <- p4 + labs(x = "Residue Index", y = "RMSF (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle = 30),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(rmsf_all[rmsf_all$variable == c("Wild-type", mutationlist[5]), ],
             aes(x=resid, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = c(1, seq(50, 500, 50), 530))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 50, 2))
p5 <- p5 + labs(x = "Residue Index", y = "RMSF (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle = 30),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/RMSF.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("Root-mean-square Fluctuation (Å)", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Residue Index", face = "bold", size = 20))
dev.off()

# ROG --------------
rog_wt_r1 <- read.table("3_WT/3_production/rog_seokWT.dat", sep = ",")
rog_wt_r2 <- read.table("3_WT/4_production_repeat2/rog_seokWT.dat", sep = ",")
rog_s50p_r1 <- read.table("7_S50P/3_production/rog_seokS50P.dat", sep = ",")
rog_s50p_r2 <- read.table("7_S50P/4_production_repeat2/rog_seokS50P.dat", sep = ",")
rog_r266q_r1 <- read.table("8_R266Q/3_production/rog_seokR266Q.dat", sep = ",")
rog_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/rog_seokR266Q.dat", sep = ",")
rog_g317r_r1 <- read.table("5_G317R/3_production/rog_seokG317R.dat", sep = ",")
rog_g317r_r2 <- read.table("5_G317R/4_production_repeat2/rog_seokG317R.dat", sep = ",")
rog_p389l_r1 <- read.table("4_P389L/3_production/rog_seokP389L.dat", sep = ",")
rog_p389l_r2 <- read.table("4_P389L/4_production_repeat2/rog_seokP389L.dat", sep = ",")
rog_l430p_r1 <- read.table("6_L430P/3_production/rog_seokL430P.dat", sep = ",")
rog_l430p_r2 <- read.table("6_L430P/4_production_repeat2/rog_seokL430P.dat", sep = ",")

rog_r1 <- cbind(rog_wt_r1, rog_s50p_r1[2], rog_r266q_r1[2], rog_g317r_r1[2],
                rog_p389l_r1[2], rog_l430p_r1[2])
rog_r1[1] <- seq(0, 249.99, 0.1)
colnames(rog_r1) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                      "p.Pro389Leu", "p.Leu430Pro")
rog_r1 <- reshape2::melt(rog_r1[1:7], id.var = "time")
rog_r1$rep <- "Repeat 1"

rog_r2 <- cbind(rog_wt_r2, rog_s50p_r2[2], rog_r266q_r2[2], rog_g317r_r2[2],
                rog_p389l_r2[2], rog_l430p_r2[2])
rog_r2[1] <- seq(0, 249.99, 0.1)
colnames(rog_r2) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                      "p.Pro389Leu", "p.Leu430Pro")
rog_r2 <- reshape2::melt(rog_r2[1:7], id.var = "time")
rog_r2 $rep <- "Repeat 2"

rog_all <- rbind(rog_r1, rog_r2)

p1 <- ggplot(rog_all[rog_all$variable == c("Wild-type", mutationlist[1]), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 300, 25))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 50, 1))
p1 <- p1 + labs(x = "Time (ns)", y = "rog (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(rog_all[rog_all$variable == c("Wild-type", mutationlist[2]), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 300, 25))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 50, 1))
p2 <- p2 + labs(x = "Time (ns)", y = "rog (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(rog_all[rog_all$variable == c("Wild-type", mutationlist[3]), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 300, 25))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 50, 1))
p3 <- p3 + labs(x = "Time (ns)", y = "rog (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(rog_all[rog_all$variable == c("Wild-type", mutationlist[4]), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 300, 25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 50, 1))
p4 <- p4 + labs(x = "Time (ns)", y = "rog (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(rog_all[rog_all$variable == c("Wild-type", mutationlist[5]), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 300, 25))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 50, 1))
p5 <- p5 + labs(x = "Time (ns)", y = "rog (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.10, 0.9),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/Rg.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("Radius of Gyration (Å)", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 20))
dev.off()

# Stability --------------
stability_wt_r1 <- read.table("3_WT/3_production/stability_seokWT.dat", sep = ",")
stability_wt_r2 <- read.table("3_WT/4_production_repeat2/stability_seokWT.dat", sep = ",")
stability_s50p_r1 <- read.table("7_S50P/3_production/stability_seokS50P.dat", sep = ",")
stability_s50p_r2 <- read.table("7_S50P/4_production_repeat2/stability_seokS50P.dat", sep = ",")
stability_r266q_r1 <- read.table("8_R266Q/3_production/stability_seokR266Q.dat", sep = ",")
stability_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/stability_seokR266Q.dat", sep = ",")
stability_g317r_r1 <- read.table("5_G317R/3_production/stability_seokG317R.dat", sep = ",")
stability_g317r_r2 <- read.table("5_G317R/4_production_repeat2/stability_seokG317R.dat", sep = ",")
stability_p389l_r1 <- read.table("4_P389L/3_production/stability_seokP389L.dat", sep = ",")
stability_p389l_r2 <- read.table("4_P389L/4_production_repeat2/stability_seokP389L.dat", sep = ",")
stability_l430p_r1 <- read.table("6_L430P/3_production/stability_seokL430P.dat", sep = ",")
stability_l430p_r2 <- read.table("6_L430P/4_production_repeat2/stability_seokL430P.dat", sep = ",")

stability_r1 <- cbind(stability_wt_r1, stability_s50p_r1[2], stability_r266q_r1[2], stability_g317r_r1[2],
                      stability_p389l_r1[2], stability_l430p_r1[2])
# stability_r1[1] <- seq(0, 249.99, 4.901961)
stability_r1[1] <- seq(0, 250, 5)
colnames(stability_r1) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                            "p.Pro389Leu", "p.Leu430Pro")

stability_r1 <- reshape2::melt(stability_r1[1:7], id.var = "time")
stability_r1$rep <- "Repeat 1"

stability_r2 <- cbind(stability_wt_r2, stability_s50p_r2[2], stability_r266q_r2[2], stability_g317r_r2[2],
                      stability_p389l_r2[2], stability_l430p_r2[2])
stability_r2[1] <- seq(0, 250, 5)
colnames(stability_r2) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                            "p.Pro389Leu", "p.Leu430Pro")

stability_r2 <- reshape2::melt(stability_r2[1:7], id.var = "time")
stability_r2$rep <- "Repeat 2"

stability_all <- rbind(stability_r1, stability_r2)

p1 <- ggplot(stability_all[stability_all$variable %in% c("Wild-type", mutationlist[1]), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 250, 25))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 2500, 20))
p1 <- p1 + labs(x = "Time (ns)", y = "stability (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(stability_all[stability_all$variable %in% c("Wild-type", mutationlist[2]), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 250, 25))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 2500, 20))
p2 <- p2 + labs(x = "Time (ns)", y = "stability (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(stability_all[stability_all$variable %in% c("Wild-type", mutationlist[3]), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 250, 25))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 2500, 20))
p3 <- p3 + labs(x = "Time (ns)", y = "stability (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(stability_all[stability_all$variable %in% c("Wild-type", mutationlist[4]), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 250, 25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 2500, 20))
p4 <- p4 + labs(x = "Time (ns)", y = "stability (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(stability_all[stability_all$variable %in% c("Wild-type", mutationlist[5]), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 250, 25))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 2500, 20))
p5 <- p5 + labs(x = "Time (ns)", y = "stability (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/Stability.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("FoldX Stability (kcal/mol)", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 20))
dev.off()

# Stability boxplot --------------
stability_wt_r1 <- read.table("3_WT/3_production/stability_seokWT.dat", sep = ",")
stability_wt_r2 <- read.table("3_WT/4_production_repeat2/stability_seokWT.dat", sep = ",")
stability_s50p_r1 <- read.table("7_S50P/3_production/stability_seokS50P.dat", sep = ",")
stability_s50p_r2 <- read.table("7_S50P/4_production_repeat2/stability_seokS50P.dat", sep = ",")
stability_r266q_r1 <- read.table("8_R266Q/3_production/stability_seokR266Q.dat", sep = ",")
stability_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/stability_seokR266Q.dat", sep = ",")
stability_g317r_r1 <- read.table("5_G317R/3_production/stability_seokG317R.dat", sep = ",")
stability_g317r_r2 <- read.table("5_G317R/4_production_repeat2/stability_seokG317R.dat", sep = ",")
stability_p389l_r1 <- read.table("4_P389L/3_production/stability_seokP389L.dat", sep = ",")
stability_p389l_r2 <- read.table("4_P389L/4_production_repeat2/stability_seokP389L.dat", sep = ",")
stability_l430p_r1 <- read.table("6_L430P/3_production/stability_seokL430P.dat", sep = ",")
stability_l430p_r2 <- read.table("6_L430P/4_production_repeat2/stability_seokL430P.dat", sep = ",")

stability_r1 <- cbind(stability_wt_r1, stability_s50p_r1[2], stability_r266q_r1[2], stability_g317r_r1[2],
                      stability_p389l_r1[2], stability_l430p_r1[2])
stability_r1[1] <- seq(0, 249.99, 4.901961)
colnames(stability_r1) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                            "p.Pro389Leu", "p.Leu430Pro")

stability_r1 <- reshape2::melt(stability_r1[1:7], id.var = "time")
stability_r1$rep <- "Repeat 1"

stability_r2 <- cbind(stability_wt_r2, stability_s50p_r2[2], stability_r266q_r2[2], stability_g317r_r2[2],
                      stability_p389l_r2[2], stability_l430p_r2[2])
stability_r2[1] <- seq(0, 249.99, 4.901961)
colnames(stability_r2) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                            "p.Pro389Leu", "p.Leu430Pro")

stability_r2 <- reshape2::melt(stability_r2[1:7], id.var = "time")
stability_r2 $rep <- "Repeat 2"

stability_all <- rbind(stability_r1, stability_r2)

p <- ggplot(stability_all, aes(x = variable, y=value, colour=variable)) 
p <- p + geom_boxplot()
p <- p + geom_point()
p <- p + scale_y_continuous(breaks = seq(-100, 800, 10))
p <- p + labs(x = "Variant", y = "FoldX Stability (kcal/mol)")
p <- p + ggsci::scale_color_jama()
p <- p + facet_grid(~ rep)
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_blank(),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5, hjust = 1),
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = c(.85, 1.01),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p <- p + guides(color = guide_legend(override.aes = list(size = .8)))
p

pdf("_outputs/Stability_box.pdf", width = 12, height = 6)
p
dev.off()

# SASA --------------
SASA_wt_r1 <- read.table("3_WT/3_production/SASA_seokWT.dat", sep = ",")
SASA_wt_r2 <- read.table("3_WT/4_production_repeat2/SASA_seokWT.dat", sep = ",")
SASA_s50p_r1 <- read.table("7_S50P/3_production/SASA_seokS50P.dat", sep = ",")
SASA_s50p_r2 <- read.table("7_S50P/4_production_repeat2/SASA_seokS50P.dat", sep = ",")
SASA_r266q_r1 <- read.table("8_R266Q/3_production/SASA_seokR266Q.dat", sep = ",")
SASA_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/SASA_seokR266Q.dat", sep = ",")
SASA_g317r_r1 <- read.table("5_G317R/3_production/SASA_seokG317R.dat", sep = ",")
SASA_g317r_r2 <- read.table("5_G317R/4_production_repeat2/SASA_seokG317R.dat", sep = ",")
SASA_p389l_r1 <- read.table("4_P389L/3_production/SASA_seokP389L.dat", sep = ",")
SASA_p389l_r2 <- read.table("4_P389L/4_production_repeat2/SASA_seokP389L.dat", sep = ",")
SASA_l430p_r1 <- read.table("6_L430P/3_production/SASA_seokL430P.dat", sep = ",")
SASA_l430p_r2 <- read.table("6_L430P/4_production_repeat2/SASA_seokL430P.dat", sep = ",")

SASA_r1 <- cbind(SASA_wt_r1, SASA_s50p_r1[2], SASA_r266q_r1[2], SASA_g317r_r1[2],
                 SASA_p389l_r1[2], SASA_l430p_r1[2])
SASA_r1[1] <- seq(0, 249.99, 0.1)
colnames(SASA_r1) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
SASA_r1 <- reshape2::melt(SASA_r1[1:7], id.var = "time")
SASA_r1$rep <- "Repeat 1"

SASA_r2 <- cbind(SASA_wt_r2, SASA_s50p_r2[2], SASA_r266q_r2[2], SASA_g317r_r2[2],
                 SASA_p389l_r2[2], SASA_l430p_r2[2])
SASA_r2[1] <- seq(0, 249.99, 0.1)
colnames(SASA_r2) <- c("time","Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                       "p.Pro389Leu", "p.Leu430Pro")
SASA_r2 <- reshape2::melt(SASA_r2[1:7], id.var = "time")
SASA_r2 $rep <- "Repeat 2"

SASA_all <- rbind(SASA_r1, SASA_r2)

p1 <- ggplot(SASA_all[SASA_all$variable == c("Wild-type", mutationlist[1]), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 250, 25))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 40000, 500))
p1 <- p1 + labs(x = "Time (ns)", y = "SASA (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(SASA_all[SASA_all$variable == c("Wild-type", mutationlist[2]), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 250, 25))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 40000, 500))
p2 <- p2 + labs(x = "Time (ns)", y = "SASA (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(SASA_all[SASA_all$variable == c("Wild-type", mutationlist[3]), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 250, 25))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 40000, 500))
p3 <- p3 + labs(x = "Time (ns)", y = "SASA (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(SASA_all[SASA_all$variable == c("Wild-type", mutationlist[4]), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 250, 25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 40000, 500))
p4 <- p4 + labs(x = "Time (ns)", y = "SASA (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.78, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(SASA_all[SASA_all$variable == c("Wild-type", mutationlist[5]), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 250, 25))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 40000, 500))
p5 <- p5 + labs(x = "Time (ns)", y = "SASA (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/SASA.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("Solvent Accessible Surface Area (Å²)", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 20))
dev.off()

# hbonds --------------
hbonds_wt_r1 <- read.table("3_WT/3_production/hbonds_seokWT.dat", sep = " ")
hbonds_wt_r2 <- read.table("3_WT/4_production_repeat2/hbonds_seokWT.dat", sep = " ")
hbonds_s50p_r1 <- read.table("7_S50P/3_production/hbonds_seokS50P.dat", sep = " ")
hbonds_s50p_r2 <- read.table("7_S50P/4_production_repeat2/hbonds_seokS50P.dat", sep = " ")
hbonds_r266q_r1 <- read.table("8_R266Q/3_production/hbonds_seokR266Q.dat", sep = " ")
hbonds_r266q_r2 <- read.table("8_R266Q/4_production_repeat2/hbonds_seokR266Q.dat", sep = " ")
hbonds_g317r_r1 <- read.table("5_G317R/3_production/hbonds_seokG317R.dat", sep = " ")
hbonds_g317r_r2 <- read.table("5_G317R/4_production_repeat2/hbonds_seokG317R.dat", sep = " ")
hbonds_p389l_r1 <- read.table("4_P389L/3_production/hbonds_seokP389L.dat", sep = " ")
hbonds_p389l_r2 <- read.table("4_P389L/4_production_repeat2/hbonds_seokP389L.dat", sep = " ")
hbonds_l430p_r1 <- read.table("6_L430P/3_production/hbonds_seokL430P.dat", sep = " ")
hbonds_l430p_r2 <- read.table("6_L430P/4_production_repeat2/hbonds_seokL430P.dat", sep = " ")

hbonds_r1 <- cbind(hbonds_wt_r1, hbonds_s50p_r1[2], hbonds_r266q_r1[2], hbonds_g317r_r1[2],
                   hbonds_p389l_r1[2], hbonds_l430p_r1[2])
hbonds_r1[1] <- seq(0, 249.99, 0.1)
colnames(hbonds_r1) <- c("time", "Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                         "p.Pro389Leu", "p.Leu430Pro")
hbonds_r1 <- reshape2::melt(hbonds_r1[1:7], id.var = "time")
hbonds_r1$rep <- "Repeat 1"

hbonds_r2 <- cbind(hbonds_wt_r2, hbonds_s50p_r2[2], hbonds_r266q_r2[2], hbonds_g317r_r2[2],
                   hbonds_p389l_r2[2], hbonds_l430p_r2[2])
hbonds_r2[1] <- seq(0, 249.99, 0.1)
colnames(hbonds_r2) <- c("time", "Wild-type", "p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg",
                         "p.Pro389Leu", "p.Leu430Pro")
hbonds_r2 <- reshape2::melt(hbonds_r2[1:7], id.var = "time")
hbonds_r2 $rep <- "Repeat 2"

hbonds_all <- rbind(hbonds_r1, hbonds_r2)

p1 <- ggplot(hbonds_all[hbonds_all$variable == c("Wild-type", mutationlist[1]), ],
             aes(x=time, y=value, colour=variable)) 
p1 <- p1 + geom_point(size = 0.7) 
p1 <- p1 + geom_line(size = 0.65)
p1 <- p1 + scale_x_continuous(breaks = seq(0, 250, 25))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 40000, 10))
p1 <- p1 + labs(x = "Time (ns)", y = "hbonds (Å)")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + facet_grid(~ rep)
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p2 <- ggplot(hbonds_all[hbonds_all$variable == c("Wild-type", mutationlist[2]), ],
             aes(x=time, y=value, colour=variable)) 
p2 <- p2 + geom_point(size = 0.7) 
p2 <- p2 + geom_line(size = 0.65)
p2 <- p2 + scale_x_continuous(breaks = seq(0, 250, 25))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 40000, 10))
p2 <- p2 + labs(x = "Time (ns)", y = "hbonds (Å)")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + facet_grid(~ rep)
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p3 <- ggplot(hbonds_all[hbonds_all$variable == c("Wild-type", mutationlist[3]), ],
             aes(x=time, y=value, colour=variable)) 
p3 <- p3 + geom_point(size = 0.7) 
p3 <- p3 + geom_line(size = 0.65)
p3 <- p3 + scale_x_continuous(breaks = seq(0, 250, 25))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 40000, 10))
p3 <- p3 + labs(x = "Time (ns)", y = "hbonds (Å)")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 4)])
p3 <- p3 + facet_grid(~ rep)
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.92, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p4 <- ggplot(hbonds_all[hbonds_all$variable == c("Wild-type", mutationlist[4]), ],
             aes(x=time, y=value, colour=variable)) 
p4 <- p4 + geom_point(size = 0.7) 
p4 <- p4 + geom_line(size = 0.65)
p4 <- p4 + scale_x_continuous(breaks = seq(0, 250, 25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 40000, 10))
p4 <- p4 + labs(x = "Time (ns)", y = "hbonds (Å)")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 5)])
p4 <- p4 + facet_grid(~ rep)
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.78, 0.15),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

p5 <- ggplot(hbonds_all[hbonds_all$variable == c("Wild-type", mutationlist[5]), ],
             aes(x=time, y=value, colour=variable)) 
p5 <- p5 + geom_point(size = 0.7) 
p5 <- p5 + geom_line(size = 0.65)
p5 <- p5 + scale_x_continuous(breaks = seq(0, 250, 25))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 40000, 10))
p5 <- p5 + labs(x = "Time (ns)", y = "hbonds (Å)")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 6)])
p5 <- p5 + facet_grid(~ rep)
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 15),
                 legend.position = c(0.1, 0.92),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/hbonds.pdf", width = 12, height = 16)
annotate_figure(fig, left = text_grob("# of H-bonds", face = "bold", size = 20, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 20))
dev.off()
