from Bio import *
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def A():
    with open("seq_covid.gb","w") as fd:                                                                        #ouverture du fichier
            fd.write("")                                                                                        #mise du fichier a 0
    rec_list = ["SARS-CoV-2[organism] AND refseq[filter]","SARSr-CoV RaTG13 complete genome","MP789 MT121216"]  #On crée un liste où chaque élément est une recherche pour le portail NCBI
    for i in range(len(rec_list)):                                                                              #Une boucle pour chaque itération, un élément de recherche
        Entrez.email = "exemple@toto.fr"                                                                        #email de reference
        recherche = Entrez.esearch(db="Nucleotide", term=rec_list[i])                                           #recherche sur le portail NCBI dans la base de donnée nucléotide la recherche qui correspond à l'élément de la liste
        resultat_recherche = Entrez.read(recherche)                                                             #conservation du resultat dans une variable(dictionnaire)
        recherche.close()                                                                                       #fermeture du fichier car on en a plus besoin
        id_list = resultat_recherche["IdList"]                                                                  #recupere les identifiant de sequences dans une liste
        fic_seq = Entrez.efetch(db="Nucleotide", id=id_list[0], rettype="gb")                                   #lance une recherche sur le portail NCBI avec le premier élément de la liste
        ma_seq = SeqIO.read(fic_seq, "gb")                                                                      #on lit le resultat et le stock dans une variable
        with open("seq_covid.gb","a") as fd:                                                                    #ouverture du fichier en ajout
            fd.write(ma_seq.format("gb"))                                                                       #ecriture du resultat dans notre fichier ouvert
        fic_seq.close()                                                                                         #fermeture de la variable de recherche
        
        
def B():                                                                                                #fonction qui pour chaque tour de boucle (chaque genbank) va noter les donnes dans un fichier 
    with open("info_seq_covid.txt","w") as fd:                                                          #ouverture du fichier
            fd.write("")                                                                                #mise du fichier a 0
    with open("info_seq_covid.txt", "a") as fd:                                                         #ouverture du fichier en ajout
        ma_seq = SeqIO.parse("seq_covid.gb","genbank")                                                  #lecture des genbanks dans le fichier qu'on met dans une variable
        for seq in ma_seq:                                                                              #boucle qui itere de facon a separer les fichiers genbanks
            count_gene = 0                                                                              #initialisations de plusieurs variables afin qu'a chaque tour 
            count_cg = 0                                                                                #de boucle tout revienne a 0 pour pas avoir de valeurs erronées entre chaque
            count_tot = 0                                                                               #fichiers genbanks
            fd.write("nom de l’organisme dont provient la donnée : "+seq.annotations['source'] +'\n')   #on accede au nom de l'organisme et on le note
            fd.write("Taxonomie : "+str(seq.annotations['taxonomy'])+'\n')                              #on accede a la taxonomie et on la note
            fd.write("id de la sequence : "+seq.id+'\n')                                                #on accede a l'id de la sequence et on la note
            fd.write("date de création de la donnée : "+seq.annotations['date']+'\n')                   #on accede a la date et on la note
            for features in seq.features:                                                               #creation d'une boucle pour iterer sur chaque feature de notre genbank
                if features.type == 'gene':                                                             #on regarde si le type de feature est un gene et si c'est le cas :
                    count_gene +=1                                                                      #on ajoute 1 a notre compteur de gene
                    nom_gene = features.qualifiers['gene'][0]                                           #on recupere son nom
                    gene_start = features.location.start                                                #on recupere le codon start de de gene
                    gene_stop = features.location.end                                                   #idem pour le codon stop
                    fd.write("le gene n°"+str(count_gene)+' est appelé '+nom_gene+'. codon start = '+str(gene_start)+' et codon stop = '+str(gene_stop)+'\n')  #on ecrit le n° du gene, son nom,codon start,codon stop
                if features.type == 'CDS':                                                              #on regarde si le type de feature est CDS et si c'est le cas :
                    fd.write("id de la proteine : "+features.qualifiers['protein_id'][0]+'\n')          #on recupere l'id de la proteine et on l'ecrit
            fd.write("nombre de gène : "+str(count_gene)+'\n')                                          #on ecrit le nombre de gene
            sequence = seq.seq                                                                          #on stock dans une variable la sequence
            count_cg += sequence.count('C')                                                             #on ajoute au compteur CG le nombre de C
            count_cg += sequence.count('G')                                                             #et le nombre de G
            count_tot +=  sequence.count('A')                                                           #on ajoute au compteur total le nombre de A
            count_tot +=  sequence.count('T')                                                           #et de T
            count_tot +=  sequence.count('C')                                                           #et de C
            count_tot +=  sequence.count('G')                                                           #et de G
            pourcentage = (count_cg/count_tot)*100                                                      #on calcule le pourcentage de CG dans la sequence grace a ces variables
            fd.write("pourcentage de CG : "+str(pourcentage)+'%'+'\n')                                  #et on l'ecrit
            fd.write('\n')                                                                              #on saute une ligne pour que notre fichier soit plus lisible et séparer les genbank
            
A()
B()