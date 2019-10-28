##### ZOOM2PIXEL
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/25/2019
# Establishing relationship between zoom level and pixel counts
# For Cannon EOS Rebel T6

out_path = 'figs/multipin_pilot/';

z2p = readxl::read_xlsx('figs/multipin_pilot/z2p.xlsx')
z2p$ratio <- z2p$width/z2p$height
z2p$height <- z2p$width/median(z2p$ratio)
z2p$pixcount <- z2p$width * z2p$height
z2p$pixcount <- z2p$pixcount[z2p$focalLength == 55]/z2p$pixcount

pix_cor = NULL
i = 1
for (fl in unique(z2p$focalLength)) {
  pix_cor$focalLength[i] <- fl
  pix_cor$correction[i] <- mean(z2p$pixcount[z2p$focalLength == fl])
  i <- i + 1
}
pix_cor <- data.frame(pix_cor)
write.table(pix_cor, file = sprintf('%spix_cor.csv',out_path),
            sep = ",", row.names = F)

z2p.model <- lm(formula = pixcount ~ focalLength + I(focalLength^2) + I(focalLength^3), data = z2p)
z2p.model
ggplot(z2p) +
  geom_point(aes(x = focalLength, y = pixcount)) +
  geom_point(data = data.frame(seq(18,55,1)),
             aes(x = seq(18,55,1), y = predict(z2p.model, data.frame(focalLength = seq(18,55,1)))))