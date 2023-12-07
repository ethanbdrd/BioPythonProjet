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
            count_tot +=  count_cg                                                                      #et de CG
            pourcentage = (count_cg/count_tot)*100                                                      #on calcule le pourcentage de CG dans la sequence grace a ces variables
            fd.write("pourcentage de CG : "+str(pourcentage)+'%'+'\n')                                  #et on l'ecrit
            fd.write('\n')                                                                              #on saute une ligne pour que notre fichier soit plus lisible et séparer les genbank



def C(user_input):                                                                                      #creer un fichier contenant les séquences protéiques de Spike pour chaque
                                                                                                        #genbank en fonction du parametre 'user_input' qui doit etre 'S' pour la fonction C ou 'S', 'M' ou 'N' pour 'G'
    with open("spike.fasta","w") as fd:                                                                 #ouverture du fichier
            fd.write("")                                                                                #mise du fichier a 0
    with open("spike.fasta", "a") as fd:                                                                #ouverture du fichier en ajout
        ma_seq = SeqIO.parse("seq_covid.gb","genbank")                                                  #lecture des genbanks dans le fichier qu'on met dans une variable
        for seq in ma_seq:                                                                              #boucle qui itere de facon a separer les fichiers genbanks
            organism = seq.annotations['organism']                                                      #recuperation du nom de l'organisme
            for features in seq.features:                                                               #creation d'une boucle pour iterer sur chaque feature de notre genbank
                if features.type == 'CDS':                                                              #on regarde si le type de feature est un gene et si c'est le cas :
                    if features.qualifiers['gene'][0] == user_input:                                    #on regarde si c'est le gene de l'utiliateur et si c'est le cas :
                        proteine_id = features.qualifiers['protein_id'][0]                              #on recupere la protein id dans une variable
                        translation =  features.qualifiers['translation'][0]                            #on recupere la traduction dans une variable
                        product = features.qualifiers['product'][0]                                     #on recupere le product dans une variable
                        fd.write(">"+proteine_id+" "+product+" "+organism+" "+'\n'+translation+"\n")    #on note tout dans le fichier afin de creer un fichier fasta complet



from Bio.Align.Applications import MafftCommandline                            #imporation pour utilisation MAFFT
def D():                                                                       #code copié-collé de moodle pour l'alignement
    commande = MafftCommandline(input="spike.fasta")
    myStdout, myStderr = commande()
    with open("aln-spike.fasta", 'w') as w:
        w.write(myStdout)


def E():                                                                                                    #creer un fichier qui liste les positions où les sequences sont différentes
    with open("resultatComparaison_geneS.txt","w") as fd:                                                   #ouverture du fichier
            fd.write("")                                                                                    #mise du fichier a 0
    with open("resultatComparaison_geneS.txt", "a") as fd:                                                  #ouverture du fichier en ajout
        ma_seq = list(SeqIO.parse("aln-spike.fasta","fasta"))                                               #creation d'une liste d'element fasta      
        fd.write("POSITION | HOMME | CHAUVE-SOURIS | PANGOLIN"+'\n')                                        #on écrit d'abord les colonnes dans notre fichier
        seq_homme = ma_seq[0].seq                                                                           #on recupere la sequence de l'homme dans une variable
        seq_chauve_souris = ma_seq[1].seq                                                                   #idem pour la chauve-souris
        seq_pangolin = ma_seq[2].seq                                                                        #et pour le pangolin
        for i in range(len(seq_homme)):                                                                     #en creer un boucle sur i en fonction de la taille de la sequences alignées qui itère sur chaque lettre
            if (seq_homme[i]!=seq_chauve_souris[i]) or (seq_chauve_souris[i]!=seq_pangolin[i]):             #on regarde les lettre grace a l'indice i et si les lettres alignées sont différentes :
                fd.write(str(i)+' | '+seq_homme[i]+' | '+seq_chauve_souris[i]+' | '+seq_pangolin[i]+'\n')   #on les notes dans notre fichier txt avec leur position




def F():                                                                                                                #fonction pour calculer le taux de convertion
    count_cs=0                                                                                                          #initialisation d'une variable pour compter le nombre de lettre en commun  avec homme-chauve souris
    count_p=0                                                                                                           #pareil pour homme-pangolin
    ma_seq = list(SeqIO.parse("aln-spike.fasta","fasta"))                                                               #creation d'une liste d'element fasta      
    seq_homme = ma_seq[0].seq                                                                                           #on recupere la sequence de l'homme dans une variable
    seq_chauve_souris = ma_seq[1].seq                                                                                   #idem pour la chauve-souris
    seq_pangolin = ma_seq[2].seq                                                                                        #et pour le pangolin
    for i in range(len(seq_homme)):                                                                                     #en creer un boucle sur i en fonction de la taille de la sequences alignées qui itère sur chaque lettre
        if (seq_homme[i]==seq_chauve_souris[i]):                                                                        #si les lettres alignées sont les mêmes:
            count_cs+=1                                                                                                 #on ajoute 1 a count_cs
        if (seq_homme[i]==seq_pangolin[i]):                                                                             #si les lettres alignées sont les mêmes:
            count_p+=1                                                                                                  #on ajoute 1
    resultat_taux_cs = (count_cs/len(seq_homme))*100                                                                    #a la fin de la boucle on calcule le pourcentage pour chauve souris
    resultat_taux_pang = (count_p/len(seq_homme))*100                                                                   #idem pou pangolin
    print("Taux de convertion de la proteine spike de l'homme vers la chauve souris est de ",resultat_taux_cs,"%.")     #on affiche le resultat de taux de convertion pour chauve souris
    print("Taux de convertion de la proteine spike de l'homme vers le pangolin est de ",resultat_taux_pang,"%.")        #idem pour pangolin


def G():                                                                    #fonction qui demande à l'utilisateur le gene qui veut pour executer le code en fonction de sa réponse
    end = 0                                                                 #initialisation d'une variable à 0
    while end ==0:                                                          #creation d'une boucle while. Tant que l'utilisateur ne donne pas une réponse convenable, on lui repose la question
        user_input = input('Choisissez un Gene entre S,M et N...')          #phrase qui demande le gene et permet une réponse pour l'utilisateur
        if user_input == 'M' or user_input == 'N' or user_input =='S':      #si sa reponse est correcte ...
            C(user_input)                                                   #on execute les fonctions...
            D()
            E()
            F()
            end=1                                                           #on met notre variable a 1 pour sortir du while


#Question H : les proteines presentes dans les 3 coronavirus sont tres similaires car pour chaques proteines, 
    #on a un taux de de conservation superieur a 90%, le taux le plus remarquable est celui avec le gene M avec 98% et 99% de conservation.
    #Alors oui les virus se ressemblent. Celui de l'homme et de la chauve-souris semblent etre les plus proches

def I():
    with open("spike.gb","w") as fd:                                                #ouverture du fichier
        fd.write("")                                                                #mise du fichier a 0
    ma_seq = list(SeqIO.parse("spike.fasta","fasta"))                               #creation d'une liste d'element fasta      
    rec_list = [ma_seq[0].id]+[ma_seq[1].id]+[ma_seq[2].id]                         #On crée un liste où chaque élément est une recherche pour le portail NCBI en fonction du gene choisi auparavant
    for i in range(len(rec_list)):                                                  #Une boucle : pour chaque itération, un élément de recherche
        Entrez.email = "exemple@toto.fr"                                            #email de reference
        recherche = Entrez.esearch(db="Protein", term=rec_list[i])                  #recherche sur le portail NCBI dans la base de donnée Proteine la recherche qui correspond à l'élément de la liste
        resultat_recherche = Entrez.read(recherche)                                 #conservation du resultat dans une variable(dictionnaire)
        recherche.close()                                                           #fermeture du fichier car on en a plus besoin
        id_list = resultat_recherche["IdList"]                                      #recupere les identifiant de sequences dans une liste
        fic_seq = Entrez.efetch(db="Protein", id=id_list[0], rettype="gb")          #lance une recherche sur le portail NCBI avec le premier élément de la liste
        ma_seq = SeqIO.read(fic_seq, "gb")                                          #on lit le resultat et le stock dans une variable
        with open("spike.gb","a") as fd:                                            #ouverture du fichier en ajout
            fd.write(ma_seq.format("gb"))                                           #ecriture du resultat dans notre fichier ouvert
        fic_seq.close()                                                             #fermeture de la variable de recherche

#FONCTIONS A EXECUTER

#A()
#B()
#C('S')
#D()
#E()
#F()
#G()
#I()
#J()
#K()

                                                                               
def J(fichier_voulu, format_fichier, penalite_decalage=-1, score_correspondance=1, penalite_non_correspondance=-1):                             #creation d'une fonction d'alignement en fonction d'un fichier choisi
    with open("alignement.fasta","w") as fd:                                                                                                    #ouverture du fichier
            fd.write("")                                                                                                                        #mise du fichier a 0         
    ma_seq = list(SeqIO.parse(fichier_voulu,format_fichier))                                                                                    #creation d'une liste d'element fasta
    for k in range(len(ma_seq)-1):                                                                                                              #boucle pour parcourir le fichier fasta
        seq1 = ma_seq[0].seq                                                                                                                    #2 variables qui sont les sequences a aligner
        seq2 = ma_seq[k+1].seq      
        lignes, colonnes = len(seq1) + 1, len(seq2) + 1                                                                                         #Initialiser la matrice de score
        matrice_score = [[0] * colonnes for _ in range(lignes)]
        for i in range(1, lignes):                                                                                                              #Initialiser la première colonne et la première ligne avec les pénalités de décalage
            matrice_score[i][0] = matrice_score[i-1][0] + penalite_decalage
        for j in range(1, colonnes):
            matrice_score[0][j] = matrice_score[0][j-1] + penalite_decalage
        for i in range(1, lignes):                                                                                                              #Remplir la matrice de score en utilisant la programmation dynamique
            for j in range(1, colonnes):                                                                                                        #double boucle
                correspondance = matrice_score[i-1][j-1] + (score_correspondance if seq1[i-1] == seq2[j-1] else penalite_non_correspondance)
                suppression = matrice_score[i-1][j] + penalite_decalage
                insertion = matrice_score[i][j-1] + penalite_decalage
                matrice_score[i][j] = max(correspondance, suppression, insertion) 
        align1, align2 = '', ''                                                                                                                 #Effectuer la trace-back pour récupérer l'alignement
        i, j = lignes - 1, colonnes - 1
        while i > 0 or j > 0:
            if i > 0 and matrice_score[i][j] == matrice_score[i-1][j] + penalite_decalage:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            elif j > 0 and matrice_score[i][j] == matrice_score[i][j-1] + penalite_decalage:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
            else:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
        if k==0:
            with open("alignement.fasta","a") as fd:                                                                                            #ouverture du fichier
                fd.write(f'>{ma_seq[0].description}\n')
                fd.write(align1 + '\n')
                fd.write(f'>{ma_seq[k].description}\n')
                fd.write(align2 + '\n')
        else :
            with open("alignement.fasta","a") as fd:                                                                                            #ouverture du fichier
                fd.write(f'>{ma_seq[k].description}\n')
                fd.write(align2 + '\n')
        

J('spike.fasta','fasta')

            
