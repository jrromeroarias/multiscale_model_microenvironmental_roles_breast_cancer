#Cargar la librer??a 
library(BoolNet)
#Establecer la carpeta de trabajo
#setwd()
#Cargar el script
#source("Atractores WT KO OE.r") 
#Atractores Grandes WT en tabla
write("\\documentclass[10pt,a4paper,notitlepage]{report}\n\\usepackage{tabularx}\n\\usepackage{colortbl}\n\\begin{document}\n\\section*{Wild-Type}", file="atractoresWtKoOe.tex")
atpinpol <- loadNetwork("CancerBooleanN.txt")
attrs <- getAttractors(atpinpol)
par(mfrow=c(1,length (table(sapply(attrs$attractors, function(attractor){length(attractor$involvedStates)})))))
plotAttractors(attrs)
attractorsToLaTeX(attrs, onColor="[gray]{1}", offColor="[gray]{1}", file="atractores.tex")
atrcWT <- readLines("atractores.tex",encoding="UTF-8")
write(atrcWT, file="atractoresWtKoOe.tex",append=TRUE)


#Atractores Grandes KO y OE en tabla de latex
for(g in atpinpol$genes){
	for(x in 0:1){
		if(x==0){
			write(paste("\\newpage\n\\section*{",g," KO}\n",sep=''), file="atractoresWtKoOe.tex",append=TRUE)
			}
		else{
			write(paste("\\newpage\n\\section*{",g," OE}\n",sep=''), file="atractoresWtKoOe.tex",append=TRUE)
			}
		mut <- fixGenes(atpinpol, g, x)
		attrs <- getAttractors(mut)
		attractorsToLaTeX(attrs, onColor="[gray]{1}", offColor="[gray]{1}", file="atractores.tex")
		atrcWT <- readLines("atractores.tex",encoding="UTF-8")
		write(atrcWT, file="atractoresWtKoOe.tex",append=TRUE)
		}
	}

#Para terminar el .tex
write("\n\\end{document}", file="atractoresWtKoOe.tex",append=TRUE)
atrcOK <- readLines("atractoresWtKoOe.tex",encoding="UTF-8")
atrcOK <- gsub(pattern="\\begin{table}\\[ht\\]", replacement="", x=atrcOK, perl=TRUE)
atrcOK <- gsub(pattern="\\end{table}", replacement="",x=atrcOK, perl = TRUE)
write(atrcOK, file="atractoresWtKoOe.tex")

