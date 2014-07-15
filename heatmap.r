# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
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
#mydist = function(v) {return(as.dist(100 - v))}
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
pdf(args[2])

maxi=as.numeric(args[3])
trueMax=as.numeric(args[4])
n=100 # number of steps between 2 colors
breaks <- seq(from = 0, to = maxi, length = 13)
cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))
palette=colorRampPalette(c("green", "yellow", "red", "brown", "black"))(n = 5*n-1)
breaks=c(seq(0,maxi/4,length=n), # for green
               seq(maxi/4,maxi/2,length=n), # for yellow
               seq(maxi/2,3*maxi/4,length=n), # for red
               seq(3*maxi/4,maxi,length=n), # for brown
               seq(maxi,trueMax,length=n)) # for black

#layouts for legend
mat <- matrix(c(1,0,0,2), 2) 
layout(mat, c(4,10), c(4,10))

breaksToMaxi=breaks[1:(4*n)] # prend que les breaks <=maxi
black.width=maxi/10
black.space=maxi/10


par(xpd=T,cex=1.8,mar=c(1,1,1,1))
plot(c(0,maxi+black.width+black.space),c(0,2),type="n",yaxt="n",ylab="",xlab=NULL,xaxt="n",xaxs = "i", yaxs = "i")
rect(breaksToMaxi[-length(breaksToMaxi)],0,breaksToMaxi[-1],2,col=palette,border=NA)
rect(maxi+black.space,0,maxi+black.space+black.width,2,col="black",border=NA)

ti=pretty(1:maxi)
ti=ti[ti<trueMax]
axis(1,at=c(ti,maxi+black.space+black.width/2),label=c(ti,trueMax))

# pour faire un break
rect(maxi,-0.1,maxi+black.space,2.1,col="white",border=NA)

 heatmap.2(cr3,
 trace = "none",
 dendrogram = "none",
  key = FALSE,
  Rowv = "NA",
  Colv = "NA",
 col=palette,
 breaks = breaks,
 margins=c(10,10),
 main = args[5])

# par(fig=c(0,0.3,0.7,1), new=TRUE)
 
 


