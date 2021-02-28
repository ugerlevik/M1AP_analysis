##################################################
## Project: M1AP mutations
## Script purpose: visualization of analyses
## Date: November 2, 2020
## Author: Umut Gerlevik
##################################################

# Read libraries -----------
library(ggplot2)
library(ggpubr)

# Common variables ----------
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")
palet <- ggsci::pal_jama()
colorPalet <- palet(6)

# Local H-bonds numbers -------------
hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 50_wt.dat", sep = " ")
hbonds_s50p <- read.table("_outputs/hbonds/hbonds10ofresid 50_S50P.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_s50p[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[1])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p1 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p1 <- p1 + geom_boxplot(size = 1.1, alpha = 0.75)
p1 <- p1 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p1 <- p1 + labs(x = "Time (ns)", y = "Number of H-bonds")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 12),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1

hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 50_wt_r2.dat", sep = " ")
hbonds_s50p <- read.table("_outputs/hbonds/hbonds10ofresid 50_S50P_r2.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_s50p[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[1])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p2 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p2 <- p2 + geom_boxplot(size = 1.1, alpha = 0.75)
p2 <- p2 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p2 <- p2 + labs(x = "Time (ns)", y = "Number of H-bonds")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 2)])
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2

hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 266_wt.dat", sep = " ")
hbonds_r266q <- read.table("_outputs/hbonds/hbonds10ofresid 266_R266Q.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_r266q[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[2])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p3 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p3 <- p3 + geom_boxplot(size = 1.1, alpha = 0.75)
p3 <- p3 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p3 <- p3 + labs(x = "Time (ns)", y = "Number of H-bonds")
p3 <- p3 + scale_color_manual(values = colorPalet[c(1, 3)])
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 12),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 266_wt_r2.dat", sep = " ")
hbonds_r266q <- read.table("_outputs/hbonds/hbonds10ofresid 266_R266Q_r2.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_r266q[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[2])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p4 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p4 <- p4 + geom_boxplot(size = 1.1, alpha = 0.75)
p4 <- p4 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p4 <- p4 + labs(x = "Time (ns)", y = "Number of H-bonds")
p4 <- p4 + scale_color_manual(values = colorPalet[c(1, 3)])
p4 <- p4 + theme_minimal()
p4 <- p4 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p4


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 317_wt.dat", sep = " ")
hbonds_g317r <- read.table("_outputs/hbonds/hbonds10ofresid 317_G317R.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_g317r[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[3])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p5 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p5 <- p5 + geom_boxplot(size = 1.1, alpha = 0.75)
p5 <- p5 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p5 <- p5 + labs(x = "Time (ns)", y = "Number of H-bonds")
p5 <- p5 + scale_color_manual(values = colorPalet[c(1, 4)])
p5 <- p5 + theme_minimal()
p5 <- p5 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 12),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p5


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 317_wt_r2.dat", sep = " ")
hbonds_g317r <- read.table("_outputs/hbonds/hbonds10ofresid 317_G317R_r2.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_g317r[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[3])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p6 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p6 <- p6 + geom_boxplot(size = 1.1, alpha = 0.75)
p6 <- p6 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p6 <- p6 + labs(x = "Time (ns)", y = "Number of H-bonds")
p6 <- p6 + scale_color_manual(values = colorPalet[c(1, 4)])
p6 <- p6 + theme_minimal()
p6 <- p6 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p6

hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 389_wt.dat", sep = " ")
hbonds_p389l <- read.table("_outputs/hbonds/hbonds10ofresid 389_P389L.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_p389l[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[4])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p7 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p7 <- p7 + geom_boxplot(size = 1.1, alpha = 0.75)
p7 <- p7 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p7 <- p7 + labs(x = "Time (ns)", y = "Number of H-bonds")
p7 <- p7 + scale_color_manual(values = colorPalet[c(1, 5)])
p7 <- p7 + theme_minimal()
p7 <- p7 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 12),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p7


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 389_wt_r2.dat", sep = " ")
hbonds_p389l <- read.table("_outputs/hbonds/hbonds10ofresid 389_P389L_r2.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_p389l[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[4])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p8 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p8 <- p8 + geom_boxplot(size = 1.1, alpha = 0.75)
p8 <- p8 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p8 <- p8 + labs(x = "Time (ns)", y = "Number of H-bonds")
p8 <- p8 + scale_color_manual(values = colorPalet[c(1, 5)])
p8 <- p8 + theme_minimal()
p8 <- p8 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p8


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 430_wt.dat", sep = " ")
hbonds_p389l <- read.table("_outputs/hbonds/hbonds10ofresid 430_L430P.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_p389l[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[5])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p9 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p9 <- p9 + geom_boxplot(size = 1.1, alpha = 0.75)
p9 <- p9 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p9 <- p9 + labs(x = "Time (ns)", y = "Number of H-bonds")
p9 <- p9 + scale_color_manual(values = colorPalet[c(1, 6)])
p9 <- p9 + theme_minimal()
p9 <- p9 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 12),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p9


hbonds_wt <- read.table("_outputs/hbonds/hbonds10ofresid 430_wt_r2.dat", sep = " ")
hbonds_l430p <- read.table("_outputs/hbonds/hbonds10ofresid 430_L430P_r2.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_l430p[2])
hbonds_all[1] <- seq(0, 249.99, 0.05)
colnames(hbonds_all) <-  c("time","WT",mutationlist[5])
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p10 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p10 <- p10 + geom_boxplot(size = 1.1, alpha = 0.75)
p10 <- p10 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,20))
p10 <- p10 + labs(x = "Time (ns)", y = "Number of H-bonds")
p10 <- p10 + scale_color_manual(values = colorPalet[c(1, 6)])
p10 <- p10 + theme_minimal()
p10 <- p10 + theme(axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p10

fig <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, widths = c(2,2,2,2,2,2,2,2,2,2),
                 nrow = 5, ncol = 2)

pdf("_outputs/hbondsWithin10_allCases.pdf", width = 4.5, height = 13)
annotate_figure(fig, left = text_grob("Number of H-bonds", size = 18, rot = 90, face = "bold"))
dev.off()

