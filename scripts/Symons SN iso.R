##### Bulk stable isotope analyses (SIA)
library("ggplot2")
iso<-read.csv("data/Symons iso data/Symons_Yos_bulkSIA.csv")

str(iso)

# elevation
d13C.elevation<-ggplot(iso, aes(x=d13C, y=elevation.m, color=type)) +
  scale_color_manual(values=c("orchid", "lightgoldenrod3", "mediumseagreen", "royalblue2", "coral")) +
  xlab(expression(paste(delta^{13}, C, " (\u2030, V-PDB)")))+
  ylab("Elevation (m)")+
  labs(color = "")+
  coord_cartesian(xlim=c(-45, -5), ylim=c(2400, 3500)) +
  geom_point() +
  stat_ellipse()+
  theme_classic()
ggsave("figures/Symons data plots/Bulk_13C.elevation.pdf", encod="MacRoman", height=6, width=7)

#C vs N

d13C.d15N<-ggplot(iso, aes(x=d13C, y=d15N, color=type)) +
  scale_color_manual(values=c("orchid", "lightgoldenrod3", "mediumseagreen", "royalblue2", "coral")) +
  xlab(expression(paste(delta^{13}, C, " (\u2030, V-PDB)")))+
  ylab(expression(paste(delta^{15}, N, " (\u2030, Air)")))+
  coord_cartesian(xlim=c(-45, -5), ylim=c(-6, 10)) +
  labs(color = "")+
  geom_point() +
  stat_ellipse()+
  theme_classic()
ggsave("figures/Symons data plots/Bulk_13C.15N.pdf", encod="MacRoman", height=6, width=7)
