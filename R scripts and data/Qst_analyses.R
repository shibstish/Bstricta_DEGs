#=====================================
# Qst analyses
#=====================================
# use the normalized to read counts per million dataset

library(tidyverse)

Bs_normalized <- read.csv("Bs_normalized.txt", header = T)
Bs_normalized <- read.csv("z.scores.csv")
Bs_normalized <- Bs_normalized %>% column_to_rownames("X") %>% as.matrix()
rownames(Bs_normalized) <- gsub("gene:", "", rownames(Bs_normalized))

pop_A <- as.data.frame(Bs_normalized[, c(1:3)])
pop_B <- as.data.frame(Bs_normalized[,4:6])
pop_C <- as.data.frame(Bs_normalized[, 7:9])
pop_D <- as.data.frame(Bs_normalized[, 10:12])
pop_E <- as.data.frame(Bs_normalized[, c(13,14)])

## var-between ##

#sample sizes
n_a <- 3
n_b <- 3
n_c <- 3
n_d <- 3
n_e <- 2

# individual means
a.mean <- rowMeans(pop_A)
b.mean <- rowMeans(pop_B)
c.mean <- rowMeans(pop_C)
d.mean <- rowMeans(pop_D)
e.mean <- rowMeans(pop_E)

# combining the populations into pairwise assessments
a.b <- cbind(pop_A, pop_B)
a.c <- cbind(pop_A, pop_C)
a.d <- cbind(pop_A, pop_D)
a.e <- cbind(pop_A, pop_E)
b.c <- cbind(pop_B, pop_C)
b.d <- cbind(pop_B, pop_D)
b.e <- cbind(pop_B, pop_E)
c.d <- cbind(pop_C, pop_D)
c.e <- cbind(pop_C, pop_E)
d.e <- cbind(pop_D, pop_E)

# "grand" means
a.b.mean <- rowMeans(a.b)
a.c.mean <- rowMeans(a.c)
a.d.mean <- rowMeans(a.d)
a.e.mean <- rowMeans(a.e)
b.c.mean <- rowMeans(b.c)
b.d.mean <- rowMeans(b.d)
b.e.mean <- rowMeans(b.e)
c.d.mean <- rowMeans(c.d)
c.e.mean <- rowMeans(c.e)
d.e.mean <- rowMeans(d.e)


a.b.var <- (((n_a/(n_a-1))*(a.mean - a.b.mean)^2)+(n_b/(n_b-1))*(b.mean - a.b.mean)^2)
a.c.var <- (((n_a/(n_a-1))*(a.mean - a.c.mean)^2)+(n_c/(n_c-1))*(c.mean - a.c.mean)^2)
a.d.var <- (((n_a/(n_a-1))*(a.mean - a.d.mean)^2)+(n_d/(n_d-1))*(d.mean - a.d.mean)^2)
a.e.var <- (((n_a/(n_a-1))*(a.mean - a.e.mean)^2)+(n_e/(n_e-1))*(e.mean - a.e.mean)^2)
b.c.var <- (((n_b/(n_b-1))*(b.mean - b.c.mean)^2)+(n_c/(n_c-1))*(c.mean - b.c.mean)^2)
b.d.var <- (((n_b/(n_b-1))*(b.mean - b.d.mean)^2)+(n_d/(n_d-1))*(d.mean - b.d.mean)^2)
b.e.var <- (((n_b/(n_b-1))*(b.mean - b.e.mean)^2)+(n_e/(n_e-1))*(e.mean - b.e.mean)^2)
c.d.var <- (((n_c/(n_c-1))*(c.mean - c.d.mean)^2)+(n_d/(n_d-1))*(d.mean - c.d.mean)^2)
c.e.var <- (((n_c/(n_c-1))*(c.mean - c.e.mean)^2)+(n_e/(n_e-1))*(e.mean - c.e.mean)^2)
d.e.var <- (((n_d/(n_d-1))*(d.mean - d.e.mean)^2)+(n_e/(n_e-1))*(e.mean - d.e.mean)^2)

## var-within ##

# variance within (singles)

A.var <- apply(pop_A,1,var)
B.var <- apply(pop_B,1,var)
C.var <- apply(pop_C,1,var)
D.var <- apply(pop_D,1,var)
E.var <- apply(pop_E,1,var)


### variance "withins" (pairwise) 

a.b.w <- (A.var + B.var)/6
a.c.w <- (A.var + C.var)/6
a.d.w <- (A.var + D.var)/6
a.e.w <- (A.var + E.var)/5
b.c.w <- (B.var + C.var)/6
b.d.w <- (B.var + D.var)/6
b.e.w <- (B.var + E.var)/5
c.d.w <- (C.var + D.var)/6
c.e.w <- (C.var + E.var)/5
d.e.w <- (D.var + E.var)/5

# pairwise Qst (appending qst onto the paired datasets)
a.b$qst <- (a.b.var/(a.b.var + a.b.w))
a.c$qst <- (a.c.var/(a.c.var + a.c.w))
a.d$qst <- (a.d.var/(a.d.var + a.d.w))
a.e$qst <- (a.e.var/(a.e.var + a.e.w))

b.c$qst <- (b.c.var/(b.c.var + b.c.w))
b.d$qst <- (b.d.var/(b.d.var + b.d.w))
b.e$qst <- (b.e.var/(b.e.var + b.e.w))

c.d$qst <- (c.d.var/(c.d.var + c.d.w))
c.e$qst <- (c.e.var/(c.e.var + c.e.w))

d.e$qst <- (d.e.var/(d.e.var + d.e.w))

# average genome-wide qst
gw_q <- data.frame("a-b" = a.b$qst, "a-c" = a.c$qst, "a-d" = a.d$qst, "a-e" = a.e$qst, "b-c" = b.c$qst, "b-d" = b.d$qst, "b-e" = b.e$qst, "c-d" = c.d$qst, "c-e" = c.e$qst, "d-e" = d.e$qst)
rownames(gw_q) <- rownames(a.b)


gw_q$avg <- rowMeans(gw_q)
hist(gw_q$avg) #distribution

# average qst
qst_table <- c(mean(a.b$qst, na.rm = T), mean(a.c$qst, na.rm = T), mean(a.d$qst, na.rm = T), mean(a.e$qst, na.rm = T), mean(b.c$qst, na.rm = T), mean(b.d$qst, na.rm = T), mean(b.e$qst, na.rm = T), mean(c.d$qst, na.rm = T), mean(c.e$qst, na.rm = T), mean(d.e$qst, na.rm = T))

qst_df <- as.data.frame(qst_table)
colnames(qst_df) <- "Qst"
qst_df$pop.comp <- c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342")

qst_var <- c(var(a.b$qst, na.rm = T), var(a.c$qst, na.rm = T), var(a.d$qst, na.rm = T), var(a.e$qst, na.rm = T), var(b.c$qst, na.rm = T), var(b.d$qst, na.rm = T), var(b.e$qst, na.rm = T), var(c.d$qst, na.rm = T), var(c.e$qst, na.rm = T), var(d.e$qst, na.rm = T))

qst_df$var <- qst_var



# standard errors - qst
#ab.se <- (sqrt(qst_df$Qst[1])/sqrt(6))
#ac.se <- (sqrt(qst_df$Qst[2])/sqrt(6))
#ad.se <- (sqrt(qst_df$Qst[3])/sqrt(6))
#ae.se <- (sqrt(qst_df$Qst[4])/sqrt(5))
#bc.se <- (sqrt(qst_df$Qst[5])/sqrt(6))
#bd.se <- (sqrt(qst_df$Qst[6])/sqrt(6))
#be.se <- (sqrt(qst_df$Qst[7])/sqrt(5))
#cd.se <- (sqrt(qst_df$Qst[8])/sqrt(6))
#ce.se <- (sqrt(qst_df$Qst[9])/sqrt(5))
#de.se <- (sqrt(qst_df$Qst[10])/sqrt(5))

# full data set
ab.se <- (sqrt(qst_df$Qst[1])/sqrt(nrow(qst_df)))
ac.se <- (sqrt(qst_df$Qst[2])/sqrt(nrow(qst_df)))
ad.se <- (sqrt(qst_df$Qst[3])/sqrt(nrow(qst_df)))
ae.se <- (sqrt(qst_df$Qst[4])/sqrt(nrow(qst_df)))
bc.se <- (sqrt(qst_df$Qst[5])/sqrt(nrow(qst_df)))
bd.se <- (sqrt(qst_df$Qst[6])/sqrt(nrow(qst_df)))
be.se <- (sqrt(qst_df$Qst[7])/sqrt(nrow(qst_df)))
cd.se <- (sqrt(qst_df$Qst[8])/sqrt(nrow(qst_df)))
ce.se <- (sqrt(qst_df$Qst[9])/sqrt(nrow(qst_df)))
de.se <- (sqrt(qst_df$Qst[10])/sqrt(nrow(qst_df)))

Qse.table <- data.frame("SE" = c(ab.se, ac.se, ad.se, ae.se, bc.se, bd.se, be.se, cd.se, ce.se, de.se), "pop.comp" = c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342"))

ggplot(qst_df, aes(x = pop.comp, y = Qst)) + 
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.025)+
  geom_bar(stat = "identity", position = "dodge", fill = "cadetblue")+
  geom_errorbar(aes(ymin = Qse.table$SE, ymax = Qst + Qse.table$SE), width = 0.1, color = "black")+
  labs(x="pairwise comparison", y = "Qst value")+
  coord_cartesian(ylim = c(0,1)) +
  #coord_flip(ylim = c(0,1))+
  theme_linedraw() +
  theme(legend.position="none")

#=========
# Fst 
#=========
# Fst calculated from the SNPs from RNA (Busch group)
fst_mx <- read.csv("hierFstat_pairwise_fst.csv")

fst_df <- data.frame(Fst = c(fst_mx[2,2], fst_mx[3,2],fst_mx[4,2],fst_mx[5,2],fst_mx[3,3],
                            fst_mx[4,3],fst_mx[5,3],fst_mx[4,4],fst_mx[5,4],fst_mx[5,5]))

fst_df$pop_comp <- c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342")


# SNPs-from-RNA Fst standard error
ab.Fse <- (sqrt(fst_df$Fst[1])/sqrt(6))
ac.Fse <- (sqrt(fst_df$Fst[2])/sqrt(6))
ad.Fse <- (sqrt(fst_df$Fst[3])/sqrt(6))
ae.Fse <- (sqrt(fst_df$Fst[4])/sqrt(5))
bc.Fse <- (sqrt(fst_df$Fst[5])/sqrt(6))
bd.Fse <- (sqrt(fst_df$Fst[6])/sqrt(6))
be.Fse <- (sqrt(fst_df$Fst[7])/sqrt(5))
cd.Fse <- (sqrt(fst_df$Fst[8])/sqrt(6))
ce.Fse <- (sqrt(fst_df$Fst[9])/sqrt(5))
de.Fse <- (sqrt(fst_df$Fst[10])/sqrt(5))



# for the SNPs-from-RNA Fst comparison
se.table <- data.frame("qst_se" = c(ab.se, ac.se, ad.se, ae.se, bc.se, bd.se, be.se, cd.se, ce.se, de.se),
"fst_se" = c(ab.Fse, ac.Fse, ad.Fse, ae.Fse, bc.Fse, bd.Fse, be.Fse, cd.Fse, ce.Fse, de.Fse),
"pop.comp" = c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342"))

# Fst calculated from low-coverage DNA data in ANGSD (Anderson group)

fst_df <- read.csv("newFST.csv")
fst_df[1] <- c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342")


# for the ANGSD Fst-Qst comparison
var.table <- data.frame("qst_var" = qst_var, "fst_var" = fst_df$Variance, "pop.comp" = c("2553_2710", "2553_2890", "2553_3133", "2553_3342", "2710_2890", "2710_3133", "2710_3342", "2890_3133", "2890_3342", "3133_3342"))

library(reshape2)

qst_fst <- cbind(qst_df[,c(2,1)], fst_df$Fst)
colnames(qst_fst) <- c("pop.comp", "Qst", "Fst")
qst_fst2 <- melt(qst_fst, id.vars = "pop.comp")

# For the SNPs-from-RNA
se.table2 <- melt(se.table, id.vars = "pop.comp")
se.table2$value <- sqrt(se.table2$value) 

# For the ANGSD Fst-Qst comparison
var.table2 <- melt(var.table, id.vars = "pop.comp")
var.table2$value <- sqrt(var.table2$value) # standard deviation

library(wesanderson)

# Qst-Fst comparison barplot with SE bars
ggplot(qst_fst2, aes(x = pop.comp, y = value, fill = factor(variable))) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = value-se.table2$value, ymax = value + se.table2$value), width=0.2, position=position_dodge(0.9), color = "grey41")+
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2))+
  labs(x = "pairwise comparison", y = "Fst or Qst") + 
  coord_flip(ylim = c(0,1), expand = FALSE)+
  theme(legend.position = "right") +
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold", colour = "black"), axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        axis.title.x = element_text(size=14, face="bold", colour = "black"), legend.text = element_text(size = 14, face = "bold"))




# ANGSD Fst barplot with variance bars (SUPPLEMENT)
ggplot(fst_df, aes(x = pop_comp, y = Fst)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "aquamarine4", color = "aquamarine4") +
  geom_errorbar(aes(ymin = Fst - Variance, ymax = Fst + Variance), width=0.2, position=position_dodge(0.9), color = "black")+
  labs(x = "pairwise comparison", y = "Fst") + 
  coord_flip(ylim = c(0,1), expand = FALSE)+
  theme(legend.position = "right") +
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold", colour = "black"), axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=14, face = "bold"))

# ANGSD Fst table (SUPPLEMENT)
fst_tbl <- data.frame("Populations" = fst_df$pop_comp,
                      "Fst" = round(fst_df$Fst, digits = 4),
                      "Variance" = round(fst_df$Variance, digits =4))
library(gt)
Ftbl <- fst_tbl %>% 
  gt() %>%
  tab_header(title = html("Pairwise F<sub>ST</sub>"),
             subtitle = md("Calculated from 155 individuals in ANGSD")) %>%
  cols_label(Fst = html("F<sub>ST</sub>"))

#=====================
# significance testing
# ====================
# Fst with RNA data
# randomly select 19290 instances of Fst 
gw_f <- read.csv("hierFstat_sitewise_fst.csv")
names(gw_f) <- c("X", "average_fst")

tmp <- gw_f[!is.na(gw_f$average_fst),]

sm_f <- tmp[sample(nrow(tmp), 19290, replace = F),]

gw_q_f <- data.frame(gw_q$avg, sm_f$average_fst)

### significance: mann-whitney u test

wilcox.test(gw_q_f$gw_q.avg, gw_q_f$sm_f.average_fst, conf.level = 0.01) #significant! p < 2.2e-16 




# Qst-Fst ratio
q_f <- data.frame("qst" = qst_df$Qst, "fst" = fst_df$Fst)
q_f$ratio <- q_f$qst/q_f$fst



## Violin Plots
library(wesanderson)
x <- qst_fst2
x$pop.comp <- as.factor(x$pop.comp)
x$variable <- as.factor(x$variable)

## population averages
# just violins
ggplot(x, aes(x = variable, y = value, fill = variable)) + 
  geom_violin(drop = FALSE)+
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2)) +
  labs(x = "variable", y = NULL) +
  theme_classic()
# with dots
ggplot(x, aes(x = variable, y = value, fill = variable)) + 
  geom_violin(drop = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color = "grey50", fill = "grey50") +
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2)) +
  labs(x = "variable", y = NULL) +
  theme_classic()

# with boxplots
ggplot(x, aes(x = variable, y = value, fill = variable)) + 
  geom_violin(drop = FALSE)+
  geom_boxplot(width=0.1, color = "grey40") +
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size=11, face="bold", colour = "black"),
        axis.text.y = element_text(size=11, face="bold", colour = "black"),
        legend.text = element_text(size = 11, face = "bold"))

# genome-wide
y <- melt(gw_q_f, value.name = "variable", )
colnames(y) <- c("variable", "value")

ggplot(y, aes(x = variable, y = value, fill = variable)) + 
  geom_violin(drop = FALSE)+
  #geom_dotplot(binaxis='y', dotsize = 0.2,stackdir='center',position=position_dodge(1)) +
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2)) +
  labs(x = "variable", y = NULL) +
  theme_classic()

ggplot(y, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(width = 0.1)+
  #geom_dotplot(binaxis='y', dotsize = 0.2,stackdir='center',position=position_dodge(1)) +
  scale_fill_manual(name = NULL, values = wesanderson::wes_palette(name = "Moonrise2", n = 2)) +
  labs(x = "variable", y = NULL) +
  theme_classic()
