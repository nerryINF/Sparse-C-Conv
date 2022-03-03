# Sparse-C-Conv

Code C de convolution sparse - étape 1 avant impélementation sur FPGA
## contenue des fichiers
- preproc.h : Classes utiles pour le traitement et loading des listes de voxels à partir de fichier .csv
- conv3d.h : Fonctions utiles pour la convolution (anciens push contiennent les différents essais de convolutions)
- computeunit.h : classe unité de traitement et kernet
- main.cpp : initialisation des différentes classes et calcul de la convolution

## Ouput 
Fichier conv5.csv : résultat de la convolution sous forme de liste sparse : *(postion xyz de la zone impacté, résultat de la convolution)*
