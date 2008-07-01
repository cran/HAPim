`haplomax.add` <-
function(hap.trans.pere,hap.trans.mere,perf,CD,map,position,marq.hap){

# recodage des haplotypes et des genotypes
##########################################
	 hap.pop   	=	rbind(hap.trans.pere,hap.trans.mere)
	 all.marq  	=	allele.marq(hap.pop)
	 hap.trans.mere = 	recode.hap(hap.trans.mere,all.marq)
	 hap.trans.pere = 	recode.hap(hap.trans.pere,all.marq)
	 hap.pop  	=	rbind(hap.trans.pere,hap.trans.mere)

# On calcule la distance des marques par rapport à l'origine
############################################################
	 dist.marq	=	distance.marqueurs(map)

# On calcule le nombre d'intervalles entre les marques
######################################################
	 nbre.int	=	length(map) 

# calcul du nombre de marques de l'haplotype: même nbre de marques à droite et à gauche de la position de test
##############################################################################################################
	 marq.hap.left  =      marq.hap/2

# On calcule du nombre de positions de tests
############################################
	  num.pos	= 	nbre.int-marq.hap+2

#Contrôle sur les positions de tests et redimensionnement du vecteur de positions de tests en fonction des positions cohérentes
###############################################################################################################################
# 1ère vérification: les positions rentrées par l'utilisateur doivent être comprises entre la 1ère marque et la dernière marque
# Les positions sont dans position.new
         position=sort(position) #rangement des valeurs de position par ordre croissant
         position=round(position,5) #Il faut arrondir!

         if (marq.hap.left==1) {dist.marq1=dist.marq}
         if (marq.hap.left>1) {dist.marq1=dist.marq[-c(1:(marq.hap.left-1),(length(dist.marq)-marq.hap.left+2):length(dist.marq))]}

         borne.inf=round(dist.marq1[1],5) #Il faut arrondir!
         borne.sup=round(dist.marq1[length(dist.marq1)],5) #Il faut arrondir!


         diff.left=position-borne.inf 
         diff.left
         diff.right=position-borne.sup 
         diff.right

         which=(1:length(position))[(diff.left>=0)&(diff.right<=0)]
         position.new=round(sort(position[which]),5)
         position.new

         if ((length(position.new)-marq.hap.left+1)!=length(position)) stop("error in test positions",call.=FALSE) 


#initialisation des vecteurs résultats
######################################
	 Fisher		=	rep(NA,num.pos)
	 assoc.est 	= 	rep(NA,num.pos) 
	 hap.ass.est	= 	rep(NA,num.pos)
	 param.est 	= 	matrix(NA,ncol=num.pos,nrow=2)

for (i in 1:num.pos)

        {  # Début de la boucle sur les positions de tests

              #extraction des informations pour la position  i
              nbre.all.marq	=	rep(0,marq.hap)
              all.marq.int	=	list()

              for(im in 1:(marq.hap)){
                 nbre.all.marq[im]	=	length(all.marq[[(i+im-1)]])
                 all.marq.int[[im]]	=	all.marq[[(i+im-1)]]
              }

              #calcul des structures pour la position i
               res.structure	=	structure(marq.hap,nbre.all.marq)
               cor.hap		=	corresp(hap.pop[,i:(i+marq.hap-1)],res.structure)
               pi.hap		=	pi.hap.NI(res.structure,cor.hap)
               cor.pere		=	corresp(hap.trans.pere[,i:(i+marq.hap-1)],res.structure)
               cor.mere		=	corresp(hap.trans.mere[,i:(i+marq.hap-1)],res.structure)  
               hap.assoc	=	unique(c(cor.pere$assoc,cor.mere$assoc))   
               nbre.ass		=	length(hap.assoc) 


              # initialisation des vecteurs nécessaires dépendant du nombre d'associations
              F.assoc	=	rep(NA,nbre.ass)
              val.par	=	matrix(NA,nrow=nbre.ass,ncol=2)
 
              # Début de la boucle sur toutes les associations possibles:
              if(nbre.ass>1) { # début du if sur nbre.ass
                   for (j in 1:nbre.ass){  # Début de la boucle sur toutes les associations possibles
                      assoc		=	hap.assoc[j]  
                      res		=	obj.haplomax.add(perf,CD,assoc,res.structure,pi.hap,cor.pere,cor.mere)
                      sum.res	=	summary(res)
                      temp		=	sum.res[[1]]
                      F.assoc[j]	=	temp$F[1]
                      val.par[j,2]	=	res$coefficient[2]
                      val.par[j,1]	=	temp$Mean[length(temp$Mean)]
               }  # Fin de la boucle sur les associations

              i.val=which.max(F.assoc)
              if(!is.na(i.val)){
              assoc.est[i]	=	hap.assoc[i.val]
              Fisher[i]	=	F.assoc[i.val]
              param.est[,i]	=	val.par[i.val,]

	      #retrouver les allèles de l'haplotype associé
              hap.ass.est[i]	=	retrouve.all(assoc.est[i],res.structure,all.marq.int)
              } 

             }# fin du if sur nbre.assoc

         }# fin de la boucle sur les positions de tests

#On récupère les positions de test
###################################
    nbre.marq 	= 	length(dist.marq)
    pos	 	=	dist.marq[1:(nbre.marq-marq.hap+1)]

#Regroupement des résultats sous forme de tableau
#################################################
    res		=	data.frame(pos,Fisher,hap.ass.est,t(param.est)) 


# initialisation des objets allant recevoir les résultats
##########################################################

    res2	=	data.frame(c1=rep(NA,length(position.new)),c2=rep(NA,length(position.new)),c3=rep(NA,length(position.new)),c4=rep(NA,length(position.new)),c5=rep(NA,length(position.new)))

   Fisher2	=	rep(NA,length(position.new))

   variance2	=	rep(NA,length(position.new))

   effet.Q2	=	rep(NA,length(position.new))


for (i in 1:(length(pos)-1)){ # début de la boucle sur i
   for (j in 1:length(position.new)){ # début de la boucle sur j

   #cas où la position de test (position.new[j]) est comprise entre la borne sup (pos[i+1]) et la borne inf (pos[i])de l'intervalle
   #On recalcule le test de Fisher, la variance et l'effet du QTL en pondérant par rapport à la borne sup et borne inf de l'intervalle
     if ((pos[i]<position.new[j])&(position.new[j]<pos[i+1])) {
       Fisher2[j]       =      ((pos[i+1]-position.new[j])/(pos[i+1]-pos[i]))*res[i+1,2]    +   ((position.new[j]-pos[i])/(pos[i+1]-pos[i]))*res[i,2]
       variance2[j]     =      ((pos[i+1]-position.new[j])/(pos[i+1]-pos[i]))*res[i+1,4]    +   ((position.new[j]-pos[i])/(pos[i+1]-pos[i]))*res[i,4]
       effet.Q2[j]      =      ((pos[i+1]-position.new[j])/(pos[i+1]-pos[i]))*res[i+1,5]    +   ((position.new[j]-pos[i])/(pos[i+1]-pos[i]))*res[i,5]
       res2[j,]		=      cbind(position.new[j],round(Fisher2[j],4),paste(res[i,3]),round(variance2[j],4),round(effet.Q2[j],4))
      }

   #cas où la position de test (position.new[j]) est égale à la position d'une marque (pos[i])
   #Les valeurs des tests de Fisher, de la variance et de l'effet du QTL ne changent pas
     if (pos[i]==position.new[j]) {
       Fisher2[j]       =      res[i,2]
       variance2[j]     =      res[i,4]
       effet.Q2[j]      =      res[i,5]
       res2[j,]		=      cbind(position.new[j],round(Fisher2[j],4),paste(res[i,3]),round(variance2[j],4),round(effet.Q2[j],4))
      }

   #cas où la position est supérieure à la position de l'avant dernière marque
   #On reprend les valeurs des tests de Fisher, de la variance et de l'effet du QTL de l'avant dernière marque
     if (((pos[length(pos)]<position.new[j])&(dist.marq[length(dist.marq)]>=position.new[j]))||(pos[length(pos)]==position.new[j])){
       Fisher2[j]       =      res[length(pos),2]
       variance2[j]     =      res[length(pos),4]
       effet.Q2[j]      =      res[length(pos),5]
       res2[j,]		=      cbind(position.new[j],round(Fisher2[j],4),paste(res[length(pos),3]),round(variance2[j],4),round(effet.Q2[j],4))
      }
   
    }# début de la boucle sur j

 } # fin de la boucle sur i


#On nomme les colonnes du data.frame
dimnames(res2)[[2]]=c("position","F","haplotype","variance","effect.Q")

res2

}

