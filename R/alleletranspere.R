`alleletranspere` <-
function(all1.pere,all2.pere,all1.fils,all2.fils,PLA,seuil){

n.desc= length(all1.fils)
alltrans.pere=rep(NA,n.desc)

#les deux allèles du fils sont égaux
#####################################
which=(1:n.desc)[!is.na(all1.fils) & !is.na(all2.fils) & all1.fils==all2.fils]
     alltrans.pere[which]=all1.fils[which]

#les deux allèles du père sont égaux
####################################
if(all1.pere==all2.pere){
   which=(1:n.desc)[!is.na(all1.fils) &  all1.fils==all1.pere]
      alltrans.pere[which]=all1.fils[which]
   which=(1:n.desc)[!is.na(all2.fils) &  all2.fils==all1.pere]
      alltrans.pere[which]=all2.fils[which] 
}

#les deux allèles du père ne sont pas égaux
###########################################
if(all1.pere!=all2.pere){
   # c'est l' allèle du chromosome 1 qui est transmis par le père
    which=(1:n.desc)[PLA>=seuil ]
      alltrans.pere[which]=all1.pere

# c'est l' allèle du chromosome 2 qui est transmis par le père
    which=(1:n.desc)[PLA <= (1-seuil) ]
    alltrans.pere[which]=all2.pere
}

alltrans.pere

}

