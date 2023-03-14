iso<-read.csv("isotope test/d13C_isotopes.csv")

str(iso)

# elevation
plot<-ggplot(iso, aes(x=d13C, y=elevation.m, color=type)) +
  geom_point() +
  stat_ellipse()+
  theme_bw()
ggsave("isotope test/13C.elevation.pdf", height=7, width=7)

#C vs N

plot2<-ggplot(iso, aes(x=d13C, y=d15N, color=type)) +
  geom_point() +
  stat_ellipse()+
  theme_bw()
ggsave("isotope test/13C.15N.pdf", height=7, width=7)
