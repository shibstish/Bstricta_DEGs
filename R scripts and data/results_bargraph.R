library(tidyverse)

go_res <- read_excel("./Metascape_Bstricta_05302025/metascape_result.xlsx", sheet = "Enrichment")


toPlot <- as.data.frame(subset(go_res, grepl('_Summary', go_res$GroupID)))

toPlot <- toPlot %>%
  arrange(`Log(q-value)`)

toPlot <- toPlot[1:11,]

toPlot$go_names <- paste0(toPlot$Term, ": ", toPlot$Description)

nrich <- toPlot %>% ggplot(aes(x = `Log(q-value)`, y = reorder(go_names, -`Log(q-value)`))) +
  geom_bar(stat = "identity", color = NA, fill = colorspace::heat_hcl(n = 11, h = c(0, -100), l=c(74, 40), c=c(40, 80), power = 1)) +
  geom_vline(xintercept = log10(0.05), color = "black") +
  labs(y = NULL) +
  scale_y_discrete(position = "right")+
  scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
  #scale_x_continuous() +
  theme_minimal()+
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, face = "bold", color = "black", hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

