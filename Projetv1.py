from Bio import *
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def A():
    rec_list = ["SARS-CoV-2 AND refseq[filter]","SARSr-CoV RaTG13","MP789 MT121216"]     #On crée un liste où chaque élément est une recherche pour le portail NCBI
    for i in range(len(rec_list)):                                                       #Une boucle pour chaque itération, un élément de recherche
        Entrez.email = "exemple@toto.fr"                                                 #email de reference
        recherche = Entrez.esearch(db="Nucleotide", term=rec_list[i])                    #recherche sur le portail NCBI dans la base de donnée nucléotide la recherche qui correspond à l'élément de la liste
        resultat_recherche = Entrez.read(recherche)                                      #conservation du resultat dans une variable(dictionnaire)
        recherche.close()                                                                #fermeture du fichier car on en a plus besoin
        id_list = resultat_recherche["IdList"]                                           #recupere les identifiant de sequences dans une liste
        fic_seq = Entrez.efetch(db="Nucleotide", id=id_list[0], rettype="gb")            #lance une recherche sur le portail NCBI avec le premier élément de la liste
        ma_seq = SeqIO.read(fic_seq, "gb")                                               #on lit le resultat et le stock dans une variable
        with open("seq_covid.gb","a") as fd:                                             #ouverture du fichier
            fd.write(ma_seq.format("gb"))                                                #ecriture du resultat dans notre fichier ouvert
        fic_seq.close()                                                                  #fermeture de la variable de recherche
        
        
def B():
    print('a')

A()

