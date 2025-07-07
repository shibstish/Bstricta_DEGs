library(tidyverse)

df <- data.frame(pops = c("2553-2710", "2553-2890", "2253-3133", "2553-3342", "2710-2890", "2710-3133", "2710-3342", "2890-3133", "2890-3342", "3133-3342"),
                   degs = c(152, 207, 222, 734, 165, 177, 437, 140, 253, 333),
                  ele_dist =c(2710-2553, 2890-2553, 3133-2553, 3342-2553, 2890-2710, 3133-2710, 3342-2710, 3133-2890, 3342-2890, 3342-3133))

colors <- c("hotpink", "orchid3", "deeppink3", "orangered3", "darkorange2","goldenrod", "gold", "limegreen",
            "turquoise3", "skyblue3")


df %>%
  ggplot(aes(x = ele_dist, y = degs, color = ele_dist)) +
  #geom_label(label = df$pops, nudge_y = 10, nudge_x = 10) + #keep working on nudge_y
  geom_point(size = 3) +
  scale_color_gradient(low = "orangered3", high = "turquoise3", name = "Elevational\nDistance (m)", labels = c(200, 400, 600, 800))+
  labs(x = "Distance (m) between populations", y = "No. differentially expressed genes") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12), 
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title.align=0.5)

  
