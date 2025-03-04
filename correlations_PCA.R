#PC1 = environment

# correlation between modules and environment

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#google_vignette

#===============================
# first, determine the correlation between environmental variables,
# the dominant PC, and the environmental variable correlated with the dominant PC

library(factoextra)
# environmental data on a 800m scale
geno <- c(rep("Bstricta", 14))
pop <- as.factor(c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3), rep("E", 2)))
band <- as.factor(c(rep("montane",6),rep("subapline", 3), rep("alpine", 5)))
elevation <- c(rep(2553, 3), rep(2710, 3), rep(2890,3), rep(3133, 3), rep(3342, 2))
aspect <- c(rep(140.19, 3), rep(218.65, 3), rep(66.8,3), rep(26.5, 3), rep(161.5, 2))
slope <- c(rep(11.6, 3), rep(2.8, 3), rep(6.3, 3), rep(14.9, 3), rep(16.4, 2))
water <- c(rep(3, 3), rep(1, 3), rep(5, 6), rep(3, 2)) #1 = low; 3 = moderate; 5 = high
slope_solar_4km <- c(rep(16.73, 3), rep(16.44, 3), rep(15.65, 3), rep(14.6, 5)) #collected on 07/26/2024
radiation <- c(rep(23.57, 3), rep(23.27, 3), rep(21.94, 3), rep(21.9, 3), rep(20.85, 2))#collected on 10/09/2024 (sloped solar irradiance, 800m resolution; June, July, and August values averaged. MJ/m2/day)
clear_solar800m <- c(rep(31.16, 3), rep(31.00, 3), rep(30.20, 3), rep(31.26, 3), rep(30.5, 2)) #collected on 10/09/2024 (June, July, and August values averaged. MJ/m2/day)
mxvpd_4km <- c(rep(12.09, 3), rep(11.25, 3), rep(9.46, 3), rep(8.75, 5)) #collected on 07/26/2024
mxvpd <- c(rep(23.22, 3), rep(22.16, 3), rep(20.16, 3), rep(17.00, 3), rep(14.58, 2))#collected on 10/09/2024 (800m resolution)
temperature <- c(rep(14.34, 3), rep(13.07, 3), rep(12.81, 3), rep(12.16, 3), rep(8.99, 2))#collected on 10/09/2024 

col_meta <- data.frame(geno, pop, band, aspect, slope, water, radiation, mxvpd, temperature)#, elevation)
rownames(col_meta) <- c("A10", "A26", "A5", "B19", "B28", "B30", "C10", "C11", "C22", "D14", "D16", "D20", "E13", "E7")

numbers <- unlist(lapply(col_meta, is.numeric), use.names = FALSE)
x <- col_meta[,numbers]

# correlation between environmental variables
res <- cor(x)
res <- round(res, 2)

pca_x2 <- prcomp(x, scale = T)
var <- get_pca_var(pca_x2)
# scree plot
fviz_eig(pca_x2)+labs(title = NULL) +xlab("PCs")
fviz_contrib(pca_x2, choice = "var", axes = 1, top = 10) + labs(title = "Contributions of variables to PC 1")+theme(plot.title = element_text(hjust = 0.5))
fviz_contrib(pca_x2, choice = "var", axes = 2, top = 10) + labs(title = "Contributions of variables to PC 2")+theme(plot.title = element_text(hjust = 0.5))

# variables plot
fviz_pca_var(pca_x2,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T) + xlab("PCA 1") + ylab("PCA 2") +theme(plot.title = element_text(hjust = 0.5))


pca_eigenvals <- pca_x2$x

write.csv(pca_eigenvals, "pca_eigenvals.csv")
