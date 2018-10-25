# Le script d'évaluation eval_simu.xxx (R ou Python) prendra en entrée :
# un nom de répertoire d'entrée
# un nom de répertoire de sortie
# le nombre d'exécutions de la méthode à réaliser sur chaque fichier du jeu de données simulées
# (nommé n_runs ci après)
# Pour chaque fichier <identifiant_fichier_i.txt> du jeu de données simulées, un fichier
# <identifiant_fichier_i>_results.txt de n_runs mots dans {TP, FP, FN} sera généré.


# generation de données simulées via script.
# Sortie:
#     SNP et score
#     {var1, var2, var3} <score>
#     {var1, var45, var2000, var5000} <score>
#     ...
#     {var67, var340} <score>
#
# On doit pouvoir determine Faux positifs et faux negatifs afin de pouvoir determiner recall precision.
#
# Exemple de fichier <identifiant_fichier_i>_results.txt :
# TP
# FN
# TP
# FP
# ...
# FN
#
# Regle pour pattern de taille 2
# si contient pattern simulé ca va donc TP true Positive
#
# Faux positif (on a trouve un truc faux)
# Faux negatif (on a pas trouve)
#     Fichier vide
#     Pas la bonne taille (donc la inferieur a 2)
#
# pour avoir un seul resultat par fichier, si on trouve un TP tout le fichier est TP
# sinon c'est a la majorité
#
# f measure pour eviter qu il plante a cause division par 0, si y a un seul TP et le reste en FP ou FP ca passe mais si pas de TP et tt en FN alors division par 0. Si TP=0 alors pb. Dans ce cas la j en enleverai un a celui qui peut pas etre calcule et je le met a l autre qui est a 0
