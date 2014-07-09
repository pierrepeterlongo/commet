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
breaks <- seq(from = 0, to = as.numeric(args[3]), length = 13)
cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))
# print(cr3)
#heatmap.2(cr3, dendrogram='none', Colv = "Rowv", distfun=mydist, col = topo.colors(12), scale="none", tracecol=TRUE, margins=c(10,10), key = TRUE, keysize = 1.5, denscol="1", density.info=c("histogram"), lhei = c(2, 8), breaks = breaks)

my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

col_breaks = c(seq(0,as.numeric(args[3])/2,length=100), # for red
seq(as.numeric(args[3])/2,as.numeric(args[3]),length=100), # for yellow
seq(as.numeric(args[3]),as.numeric(args[4]),length=100)) # for green

 heatmap.2(cr3,
 trace = "none",
 Rowv = "NA",
 Colv = "NA",
 col=my_palette,
 breaks = col_breaks,
 margins=c(10,10),
 main = args[5])


