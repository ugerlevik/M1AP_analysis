#########################################################################
## Project: M1AP                                 
## Script purpose: Visualize salt bridges
## Date: November 3, 2020                                                   
## Author: Umut Gerlevik                                                 
#########################################################################
library(ggplot2)
library(reshape2)

localRes_50 <- read.delim("_outputs/saltbridges/localResidueList_50.txt", header = F, sep = "}")
localRes_50 <- localRes_50[,1:(ncol(localRes_50)-1)]
localRes_50 <- gsub("[{| ]", "", localRes_50)

localRes_266 <- read.delim("_outputs/saltbridges/localResidueList_266.txt", header = F, sep = "}")
localRes_266 <- localRes_266[,1:(ncol(localRes_266)-1)]
localRes_266 <- gsub("[{| ]", "", localRes_266)

localRes_317 <- read.delim("_outputs/saltbridges/localResidueList_317.txt", header = F, sep = "}")
localRes_317 <- localRes_317[,1:(ncol(localRes_317)-1)]
localRes_317 <- gsub("[{| ]", "", localRes_317)

localRes_389 <- read.delim("_outputs/saltbridges/localResidueList_389.txt", header = F, sep = "}")
localRes_389 <- localRes_389[,1:(ncol(localRes_389)-1)]
localRes_389 <- gsub("[{| ]", "", localRes_389)

localRes_430 <- read.delim("_outputs/saltbridges/localResidueList_430.txt", header = F, sep = "}")
localRes_430 <- localRes_430[,1:(ncol(localRes_430)-1)]
localRes_430 <- gsub("[{| ]", "", localRes_430)

localRes_List <- list(localRes_50, localRes_266, localRes_317, localRes_389, localRes_430)
mutationlist <-  c("p.Ser50Pro", "p.Arg266Gln", "p.Gly317Arg", "p.Pro389Leu", "p.Leu430Pro")
mutList <- c("S50P", "R266Q", "G317R", "P389L", "L430P")
res <- c(50, 266, 317, 389, 430)
subgroup <- c(7, 8, 5, 4, 6)

for (i in 1:5) {
  saltbrFileName <- c()
  for (k in 1:length(localRes_List[[i]])) {
    saltbrFileName <- c(saltbrFileName, dir("3_WT/4_production_repeat2/saltbridges_seokWT/", pattern = ".dat")[grep(paste0("\\b",localRes_List[[i]][k],"\\b"), dir("3_WT/4_production_repeat2/saltbridges_seokWT/", pattern = ".dat"))])
    saltbrFileName <- c(saltbrFileName, dir(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"), pattern = ".dat")[grep(paste0("\\b",localRes_List[[i]][k],"\\b"), dir(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"), pattern = ".dat"))])
  }
  for (j in saltbrFileName) {
    if (j %in% dir("3_WT/4_production_repeat2/saltbridges_seokWT/") &&
        j %in% dir(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"))) {
      
      wt <- read.table(paste("3_WT/4_production_repeat2/saltbridges_seokWT/", j, sep = ""), sep = " ")
      mut <- read.table(paste0(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"), j), sep = " ")
      wt_mut <- cbind(wt, mut[2])
      wt_mut[1] <- seq(0, 249.99, 0.05)
      colnames(wt_mut) <- c("Time", "Wild-type", mutationlist[i])
      wt_mut <- melt(wt_mut[1:3], id.var = "Time")
      p <- ggplot(wt_mut, aes(x = Time, y = value, colour = variable))
      p <- p + geom_point(size = 1)
      p <- p + geom_line(size = 0.8)
      p <- p + scale_x_continuous(breaks = seq(0, 250, 25))
      p <- p + scale_y_continuous(breaks = seq(0, 300, 2))
      p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
      p <- p + ggsci::scale_color_jama()
      p <- p + theme_minimal()
      p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
                     axis.title.x = element_text(size = 24),
                     axis.ticks = element_line(),
                     axis.text = element_text(size = 19),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 15),
                     legend.position = c(.01, 1.01),
                     legend.justification = c("left", "top"),
                     legend.box.just = "left",
                     legend.margin = margin(6, 6, 6, 6),
                     panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
      pdf(file = paste0("_outputs/saltbridges/_repeat2/", mutList[i], "/WT_",mutList[i],"_", sub(".dat", "", j), ".pdf"), width = 8, height = 6)
      print(p)
      dev.off()
      
    } else if ((j %in% dir("3_WT/4_production_repeat2/saltbridges_seokWT/")) &&
               !(j %in% dir(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/")))) {
      
      wt <- read.table(paste0("3_WT/4_production_repeat2/saltbridges_seokWT/", j), sep = " ")
      wt[1] <- seq(0, 249.99, 0.05)
      colnames(wt) <- c("Time", "saltbr")
      p <- ggplot(wt, aes(x = Time, y = saltbr))
      p <- p + geom_point(size = 1)
      p <- p + geom_line(size = 0.8)
      p <- p + scale_x_continuous(breaks = seq(0, 250, 25))
      p <- p + scale_y_continuous(breaks = seq(0, 300, 2))
      p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
      p <- p + ggsci::scale_color_jama()
      p <- p + theme_minimal()
      p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
                     axis.title.x = element_text(size = 24),
                     axis.ticks = element_line(),
                     axis.text = element_text(size = 19),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 15),
                     legend.position = c(.01, 1.01),
                     legend.justification = c("left", "top"),
                     legend.box.just = "left",
                     legend.margin = margin(6, 6, 6, 6),
                     panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
      pdf(file = paste0("_outputs/saltbridges/_repeat2/", mutList[i], "/WT_", sub(".dat", "", j), ".pdf"), width = 8, height = 6)
      print(p)
      dev.off()
      
    } else if (!(j %in% dir("3_WT/4_production_repeat2/saltbridges_seokWT/")) &&
               j %in% dir(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"))) {
      
      mut <- read.table(paste0(paste0(subgroup[i],"_",mutList[i],"/4_production_repeat2/saltbridges_seok",mutList[i],"/"), j), sep = " ")
      mut[1] <- seq(0, 249.99, 0.05)
      colnames(mut) <- c("Time", "saltbr")
      p <- ggplot(mut, aes(x = Time, y = saltbr))
      p <- p + geom_point(size = 1)
      p <- p + geom_line(size = 0.8)
      p <- p + scale_x_continuous(breaks = seq(0, 250, 25))
      p <- p + scale_y_continuous(breaks = seq(0, 300, 2))
      p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
      p <- p + ggsci::scale_color_jama()
      p <- p + theme_minimal()
      p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
                     axis.title.x = element_text(size = 24),
                     axis.ticks = element_line(),
                     axis.text = element_text(size = 19),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 15),
                     legend.position = c(.01, 1.01),
                     legend.justification = c("left", "top"),
                     legend.box.just = "left",
                     legend.margin = margin(6, 6, 6, 6),
                     panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
      pdf(file = paste0("_outputs/saltbridges/_repeat2/", mutList[i], "/",mutList[i],"_", sub(".dat", "", j), ".pdf"), width = 8, height = 6)
      print(p)
      dev.off()
      
    }
  }
}

