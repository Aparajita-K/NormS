Um=as.matrix(read.table("Uselect"))
cat(file="UselectStats")
for(i in 1:dim(Um)[2])
{
    sw=shapiro.test(Um[,i])
    cat(sw$p.value,"\n",file="UselectStats",append=TRUE)
}
