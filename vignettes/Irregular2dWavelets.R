## ---- echo=TRUE----------------------------------------------------------
library(devtools)

## ---- message=FALSE------------------------------------------------------
install_github(repo = "jantonelli111/Irregular2dWavelets")
library(Irregular2dWavelets)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
data(BostonPMsim)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
head(BostonPMsim)

## ---- echo=TRUE, eval=TRUE, fig.align='center', fig.width=4--------------
library(maps)
colors = colorRampPalette(c("yellow", "orange", "red"))
BostonPMsim$col = colors(100)[as.numeric(cut(BostonPMsim$pm,breaks = 100))]
map("state","massachusetts", xlim=range(BostonPMsim$long), ylim=range(BostonPMsim$lat))
points(BostonPMsim$long, BostonPMsim$lat, col=BostonPMsim$col, cex=.2)

