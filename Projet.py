from Bio import *
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def A():
    with open("seq_covid.gb","w") as fd:                                                 #ouverture du fichier
            fd.write("")                                                                 #mise du fichier a 0
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
    with open("info_seq_covid.txt","w") as fd:                                           #ouverture du fichier
            fd.write("")                                                                 #mise du fichier a 0
    with open("info_seq_covid.txt", "a") as fd:
        ma_seq = SeqIO.parse("seq_covid.gb","genbank")
        for seq in ma_seq:
            count_gene = 0
            posistion_gene = []
            fd.write("nom de l’organisme dont provient la donnée : "+seq.annotations['source'] +'\n')
            fd.write("Taxonomie : "+str(seq.annotations['taxonomy'])+'\n')
            fd.write("id de la sequence : "+seq.id+'\n')
            fd.write("date de création de la donnée : "+seq.annotations['date']+'\n')
            for features in seq.features:
                if features.type == 'gene':
                    count_gene +=1
                    nom_gene = features.qualifiers['gene']
                    gene_start = features.location.start
                    gene_stop = features.location.end
                    posistion_gene.append((nom_gene,gene_start,gene_stop))
                    protein_id = features.qualifiers['protein_id'][0]
            fd.write("nobre de gène : "+str(count_gene)+'\n')
            fd.write("nom du gene et position : "+str(posistion_gene)+'\n')
            fd.write("id de la protein : "+protein_id+'\n')
            fd.write('\n')
            
B()
