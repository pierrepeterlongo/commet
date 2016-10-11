# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
#   Claire LEMAITRE, claire.lemaitre@inria.fr    [06/07/16]
#
# This software is a computer program whose purpose is to find all the
# similar reads between sets of NGS reads. It also provide a similarity
# score between the two samples.
#
# Copyright (C) 2014  INRIA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
png(file=args[3],width=800,height=800,res=65)
n=100 # number of steps between 2 colors
cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))

## Computing mini-maxi for colour palette
mini=min(cr3[])
maxi=max(cr3[row(cr3)!=col(cr3)]) # ignoring the diagonal
trueMax=max(cr3[]) # typically the value in the diagonal = 100
q25=quantile(cr3[row(cr3)!=col(cr3)],0.25,1)
q50=quantile(cr3[row(cr3)!=col(cr3)],0.5,1)
q75=quantile(cr3[row(cr3)!=col(cr3)],0.75,1)

## We use the quantiles to ignore some outlier values in the matrix (values<mini will have colour of mini and values>maxi will have a colour between brown and grey23)
mini=max(q25-1.5*(q75-q25),0)
maxi=min(q75+1.5*(q75-q25),trueMax)


 breaks=unique(c(seq(mini,mini+diff/4,length=n), # for green
                seq(mini+diff/4,mini+diff/2,length=n), # for yellow
                seq(mini+diff/2,mini+3*diff/4,length=n), # for red
                seq(mini+3*diff/4,maxi,length=n), # for brown
                seq(maxi,trueMax,length=n))) # for black

palette=colorRampPalette(c("green", "yellow", "red", "brown", "grey23"))(n = length(breaks)-1)

if(trueMax.needed){
  breaks=c(seq(mini,maxi,length=4*n),seq(maxi+1e-5,trueMax,le=n))
  # breaks are equally distributed in the range mini-maxi (intervals can be different in the range maxi-trueMax, containing very few points)

} else {
  breaks=c(seq(mini,maxi,length=5*n))
}



# Dendrogram is obtained with the normalized matrix 
 cr3_norm = as.matrix(read.table(file=args[2], sep=";", header=TRUE, row.names=1))
 inv_cr3 = matrix(trueMax, ncol=dim(cr3_norm)[1], nrow=dim(cr3_norm)[1]) - cr3_norm
 distance    = dist(inv_cr3)
 cluster     = hclust(distance)
 dendrogram  = as.dendrogram(cluster)
 
 # Heatmap 
 par(fig=c(0.2,1,0,0.8))

 heatmap.2(cr3,
 trace = "none",
 dendrogram = "none",
 key = FALSE,
  Rowv=dendrogram,
  Colv = rev(dendrogram),
 col=palette,
 breaks = breaks,
 margins=c(10,10),
 main = args[4])

# Adding the colour scale
par(fig=c(0.05,0.4,0.8,1), new=TRUE)

if(trueMax.needed){

  diff=maxi-mini
  breaksToMaxi=breaks[1:(4*n)] # using only breaks from mini to maxi
  black.width=max(diff/9)
  black.space=max(diff/9)
  
  plot(c(mini,maxi+black.width+black.space),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxt="n",xaxs = "i", yaxs = "i")
  rect(breaksToMaxi[-length(breaksToMaxi)],0,breaksToMaxi[-1],2,col=palette,border=NA)

  # x axis and labels
  ti=pretty(breaksToMaxi)
  ti=ti[ti<maxi]
  trueMax
  axis(1,at=c(ti,maxi+black.space+black.width/2),label=c(ti,trueMax))
  
  # plotting the TrueMax colour separated from the rest of the plot by a white space
  rect(maxi+black.space,0,maxi+black.space+black.width,2,col=palette[5*n-1],border=NA)
  rect(maxi,-0.1,maxi+black.space,2.1,col="white",border=NA)

} else{
 plot(range(breaks),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxs = "i", yaxs = "i")
 rect(breaks[-length(breaks)],0,breaks[-1],2,col=palette,border=NA)
}

d=dev.off()


