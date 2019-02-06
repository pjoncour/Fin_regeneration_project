#### Fisher matrizome
## Pauline
## 04/02/19

all=read.table("for_fisher.txt", sep="\t", header=T, dec=",",quote="",comment.char = "")

# debadup=all[all$Ad_Deb=="up",]
# debaddown=all[all$Ad_Deb=="down",]
debad=all[all$Ad_Deb!=0,]
# debembup=all[all$Emb_Deb=="up",]
# debembdown=all[all$Emb_Deb=="down",]
debemb=all[all$Emb_Deb!=0,]

# finadup=all[all$Ad_Fin=="up",]
# finaddown=all[all$Ad_Fin=="down",]
finad=all[all$Ad_Fin!=0,]
# finembup=all[all$Emb_Fin=="up",]
# finembdown=all[all$Emb_Fin=="down",]
finemb=all[all$Emb_Fin!=0,]

## Test exact de Fisher sous-liste
B3=nrow(debemb) # Nombre total de genes dans la sous_liste
B=nrow(all) # Nombre total de genes
m3=nrow(debemb[debemb$Division!=0,]) # Nombre total de genes de la categorie dans la sous_liste
m=sum(all$Division!=0) # Nombre total de genes de la categorie
  
MonFisher = function(m3,m, B3,B) { 
  
  tab=matrix(c(m3,B3-m3,m-m3, B-B3-m+m3), nrow=2)
  z=fisher.test(tab, alternative="greater") 
  z=z$p.value 
  return(z) 
}

MonFisher(m3,m,B3,B)

library(gplots)

VENN=venn(list(ad_mat=debad$GeneID[debad$Division!=0],emb_mat=debemb$GeneID[debemb$Division!=0]))

table(debad$Division[debad$GeneID%in%debemb$GeneID])
table(debad$Division[!debad$GeneID%in%debemb$GeneID])
table(debad$Division[!debemb$GeneID%in%debad$GeneID])


VENN=venn(list(ad_mat=finad$GeneID[finad$Division!=0],emb_mat=finemb$GeneID[finemb$Division!=0]))

table(finad$Division[finad$GeneID%in%finemb$GeneID])
table(finad$Division[!finad$GeneID%in%finemb$GeneID])
table(finad$Division[!finemb$GeneID%in%finad$GeneID])


VENN=venn(list(deb_ad_mat=debad$GeneID[debad$Division!=0],deb_emb_mat=debemb$GeneID[debemb$Division!=0],
               fin_emb_mat=finemb$GeneID[finemb$Division!=0],fin_ad_mat=finad$GeneID[finad$Division!=0]))



common=attr(VENN, "intersection")$`deb_ad_mat:deb_emb_mat:fin_emb_mat:fin_ad_mat`

matcom=all[all$GeneID%in%common,]



## Test exact de Fisher 2 listes
A= # Nombre total de genes dans la liste A
A1= # Nombre de genes de la categorie dans la liste A
B= # Nombre total de genes dans la liste B
B1= # Nombre de genes de la categorie dans la liste B

MonFisher = function(A, A1, B, B1) { 
  
  tab=matrix(c(A1,B1,A-A1, B-B1), nrow=2)
  z=fisher.test(tab, "two.sided") 
  z=z$p.value 
  return(z) 
}

MonFisher(A,A1,B,B1)






