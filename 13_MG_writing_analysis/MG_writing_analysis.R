rm(list = ls())
source("../project_support.r")

dir_init("./temp")

############

# Load data
# ImpDatRepl.csv was downloaded from https://github.com/pesavage/moralizing-gods on 2019/09/04
ImpData <- read.csv("ImpDatRepl.csv", stringsAsFactors = FALSE)

# Find when writing first appeared for each NGA
writing <- ImpData %>%
  filter(Writing == 1) %>%
  group_by(NGA) %>%
  filter(row_number()==1) %>%
  dplyr::select(NGA, Time, Writing) %>%
  dplyr::rename(`Writing First Recorded` = Time)

# Find when moralising gods first appeared for each NGA
MG <- ImpData %>%
  filter(MoralisingGods == 1) %>%
  group_by(NGA) %>%
  filter(row_number()==1) %>%
  dplyr::select(NGA, Time, MoralisingGods) %>%
  dplyr::rename(`Moralizing Gods First Recorded` = Time)

# join writing and moralising gods data
MGwriting <- left_join(writing, MG)

# remove Valley of Oaxaca which is always NA for MGs
MGwriting <- MGwriting %>%
  filter(NGA != "Valley of Oaxaca") %>%
  # calculate time difference between writing and MGs first appearing (negative means MGs appeared first)
  mutate(TimeDifference = `Moralizing Gods First Recorded` -  `Writing First Recorded`) %>%
  # add whether NGA was used in analysis or not
  mutate(`Social complexity data available \nboth before and after the appearance \nof MG` = if_else(NGA == "Deccan" | NGA == "Kachi Plain" | NGA == "Kansai" | NGA == "Konya Plain" | NGA == "Latium" | NGA == "Middle Yellow River Valley" | NGA == "Niger Inland Delta" | NGA == "Orkhon Valley" | NGA == "Paris Basin" | NGA == "Sogdiana" | NGA == "Susiana" | NGA == "Upper Egypt", "Yes", "No"))

# extract just NGAs used in analysis
MGwritinganl <- MGwriting[MGwriting$`Social complexity data available 
                          both before and after the appearance 
                          of MG` == "Yes",]

# The first appearance of writing and moralizing gods across NGAs. 
ggplot(MGwriting, aes(x = `Writing First Recorded`, y = `Moralizing Gods First Recorded`, colour = `Social complexity data available \nboth before and after the appearance \nof MG`)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray50', size = 1) +
  geom_abline(slope = 1, intercept = 100, color = 'gray60', linetype = "dashed") +
  geom_abline(slope = 1, intercept = -100, color = 'gray60', linetype = "dashed") +
  geom_point(size = 3) +
  geom_text(aes(x = `Writing First Recorded`, y = `Moralizing Gods First Recorded`, label=ifelse(NGA == "Susiana" | NGA == "Kachi Plain", NGA,'')), vjust = -1.25, inherit.aes = FALSE, size = 4.5, color = "#0072B2") +
  scale_x_continuous(limits = c(-4000,2000), breaks=c(-4000,-2000,0, 2000), labels=c("4000 BCE", "2000 BCE", "0 CE", "2000 CE")) +
  scale_y_continuous(limits = c(-4000,2000), breaks=c(-4000,-2000,0, 2000), labels=c("4000 BCE", "2000 BCE", "0 CE", "2000 CE")) +
  labs(color = "Social complexity data available \nboth before and after the appearance \nof MG") +
  xlab("Writing First Recorded") +
  ylab("Moralizing Gods First Recorded") +
  scale_color_manual(values=c("#E69F00", "#0072B2")) + 
  theme_bw() +
  theme(
    legend.position = c(0.7, 0.2),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text = element_text(colour = "black", size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
ggsave("./temp/WritingMG.png", height = 6.5, width = 6.5)
dev.off()

############

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")