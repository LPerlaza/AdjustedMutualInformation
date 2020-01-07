#How many comparisons per core?

T=1000000 #Total numbre of snps
splits=180 #Total number of cores
comp=c(T:1)#Vector of comparisons that each snp has to do 
x=c(1:T)#Vector of the position of the snps
ncomp=sum(as.numeric(comp))#Total of comparisons
comp_core=ncomp/splits #number of comparisons per core

count=data.frame(Linea=c(1:T))
count$Comparisons=abs(count$Linea-T)

breaks=seq(1,sum(count$Comparisons),length.out=splits)
count$sumComp=cumsum(count$Comparisons)
count$intervalsStart=cut(count$sumComp,breaks=breaks,label=round(breaks[-splits]))
count$intervalsEnd=cut(count$sumComp,breaks=breaks,label=round(breaks[-1]-1))

x=ddply(count,~intervalsStart,summarize,min=min(as.numeric(as.character(Linea))),max=max(as.numeric(as.character(Linea))),NumCompariosns=sum(Comparisons))
