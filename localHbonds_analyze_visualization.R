##################################################
## Project: M1AP mutations
## Script purpose: visualization of local H-bonds
## Date: Dec 24, 2021
## Author: Umut Gerlevik
##################################################

# Read libraries -----------
library(ggplot2)
library(ggpubr)

# Common variables ----------
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")
colorPalet <- grDevices::colors()[grDevices::colors() %in% c("steelblue3", "royalblue2",
                                                             "orangered3", "indianred3")]

# Local H-bonds numbers -------------
hbonds_wt_s50p_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 50_wt.dat", sep = " ")
hbonds_s50p_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 50_S50P.dat", sep = " ")
hbonds_wt_s50p_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 50_wt_r2.dat", sep = " ")
hbonds_s50p_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 50_S50P_r2.dat", sep = " ")

hbonds_all_s50p <- data.frame(hbonds_wt_s50p_r1, hbonds_wt_s50p_r2[2], hbonds_s50p_r1[2], hbonds_s50p_r2[2])
hbonds_all_s50p[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_all_s50p) <-  c("time", "Wild-type_R1", "Wild-type_R2",
                           paste0(mutationlist[1], "_R1"),
                           paste0(mutationlist[1], "_R2"))
hbonds_all_s50p <- reshape2::melt(hbonds_all_s50p[1:5], id.var = "time") 

p1 <- ggplot(hbonds_all_s50p, aes(x=variable, y=value, colour=variable))
p1 <- p1 + geom_violin(size = 0.75, alpha = 0.75, draw_quantiles = .5)
# p1 <- p1 + scale_x_continuous(breaks = seq(0, 500, 50))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 22))
p1 <- p1 + labs(x = "Time (ns)", y = "Number of Local H-bonds")
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
                 legend.position = "none",
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))


hbonds_wt_r266q_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 266_wt.dat", sep = " ")
hbonds_r266q_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 266_R266Q.dat", sep = " ")
hbonds_wt_r266q_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 266_wt_r2.dat", sep = " ")
hbonds_r266q_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 266_R266Q_r2.dat", sep = " ")

hbonds_all_r266q <- data.frame(hbonds_wt_r266q_r1, hbonds_wt_r266q_r2[2], hbonds_r266q_r1[2], hbonds_r266q_r2[2])
hbonds_all_r266q[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_all_r266q) <-  c("time", "Wild-type_R1", "Wild-type_R2",
                                paste0(mutationlist[2], "_R1"),
                                paste0(mutationlist[2], "_R2"))
hbonds_all_r266q <- reshape2::melt(hbonds_all_r266q[1:5], id.var = "time") 

p2 <- ggplot(hbonds_all_r266q, aes(x=variable, y=value, colour=variable))
p2 <- p2 + geom_violin(size = 0.75, alpha = 0.75, draw_quantiles = .5)
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 500, 50))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 22))
p2 <- p2 + labs(x = "Time (ns)", y = "Number of Local H-bonds")
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
                 legend.position = "none",
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))


hbonds_wt_g317r_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 317_wt.dat", sep = " ")
hbonds_g317r_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 317_G317R.dat", sep = " ")
hbonds_wt_g317r_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 317_wt_r2.dat", sep = " ")
hbonds_g317r_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 317_G317R_r2.dat", sep = " ")

hbonds_all_g317r <- data.frame(hbonds_wt_g317r_r1, hbonds_wt_g317r_r2[2], hbonds_g317r_r1[2], hbonds_g317r_r2[2])
hbonds_all_g317r[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_all_g317r) <-  c("time", "Wild-type_R1", "Wild-type_R2",
                                 paste0(mutationlist[3], "_R1"),
                                 paste0(mutationlist[3], "_R2"))
hbonds_all_g317r <- reshape2::melt(hbonds_all_g317r[1:5], id.var = "time") 

p3 <- ggplot(hbonds_all_g317r, aes(x=variable, y=value, colour=variable))
p3 <- p3 + geom_violin(size = 0.75, alpha = 0.75, draw_quantiles = .5)
# p3 <- p3 + scale_x_continuous(breaks = seq(0, 500, 50))
p3 <- p3 + scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 22))
p3 <- p3 + labs(x = "Time (ns)", y = "Number of Local H-bonds")
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
                 legend.position = "none",
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))


hbonds_wt_p389l_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 389_wt.dat", sep = " ")
hbonds_p389l_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 389_P389L.dat", sep = " ")
hbonds_wt_p389l_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 389_wt_r2.dat", sep = " ")
hbonds_p389l_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 389_P389L_r2.dat", sep = " ")

hbonds_all_p389l <- data.frame(hbonds_wt_p389l_r1, hbonds_wt_p389l_r2[2], hbonds_p389l_r1[2], hbonds_p389l_r2[2])
hbonds_all_p389l[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_all_p389l) <-  c("time", "Wild-type_R1", "Wild-type_R2",
                                 paste0(mutationlist[4], "_R1"),
                                 paste0(mutationlist[4], "_R2"))
hbonds_all_p389l <- reshape2::melt(hbonds_all_p389l[1:5], id.var = "time") 

p4 <- ggplot(hbonds_all_p389l, aes(x=variable, y=value, colour=variable))
p4 <- p4 + geom_violin(size = 0.75, alpha = 0.75, draw_quantiles = .5)
# p4 <- p4 + scale_x_continuous(breaks = seq(0, 500, 50))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 22))
p4 <- p4 + labs(x = "Time (ns)", y = "Number of Local H-bonds")
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
                 legend.position = "none",
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))


hbonds_wt_l430p_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 430_wt.dat", sep = " ")
hbonds_l430p_r1 <- read.table("_outputs/hbonds/hbonds10ofresid 430_L430P.dat", sep = " ")
hbonds_wt_l430p_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 430_wt_r2.dat", sep = " ")
hbonds_l430p_r2 <- read.table("_outputs/hbonds/hbonds10ofresid 430_L430P_r2.dat", sep = " ")

hbonds_all_l430p <- data.frame(hbonds_wt_l430p_r1, hbonds_wt_l430p_r2[2], hbonds_l430p_r1[2], hbonds_l430p_r2[2])
hbonds_all_l430p[1] <- seq(0, 499.99, 0.1)
colnames(hbonds_all_l430p) <-  c("time", "Wild-type_R1", "Wild-type_R2",
                                 paste0(mutationlist[5], "_R1"),
                                 paste0(mutationlist[5], "_R2"))
hbonds_all_l430p <- reshape2::melt(hbonds_all_l430p[1:5], id.var = "time") 

p5 <- ggplot(hbonds_all_l430p, aes(x=variable, y=value, colour=variable))
p5 <- p5 + geom_violin(size = 0.75, alpha = 0.75, draw_quantiles = .5)
# p5 <- p5 + scale_x_continuous(breaks = seq(0, 500, 50))
p5 <- p5 + scale_y_continuous(breaks = seq(0, 50, 2), limits = c(0, 22))
p5 <- p5 + labs(x = "Time (ns)", y = "Number of Local H-bonds")
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
                 legend.position = "none",
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

fig <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)

pdf("_outputs/localHbonds_500ns.pdf", width = 5, height = 12)
annotate_figure(fig, left = text_grob("Number of Local H-bonds", face = "bold", size = 12, rot = 90),
                bottom = text_grob("Time (ns)", face = "bold", size = 12))
dev.off()


