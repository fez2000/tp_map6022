# tp_map6022

Prédiction du risque de départ d'un employer a l'aide du dataset [turnover.csv](./datasets/turnover.csv) et [turnover_prepared.csv](./datasets/turnover_prepared.csv). Vous trouverez si join le [code source](./TP_MAP6009.ipynb) des expérimentations du [Rapport]() .

- [Description du dataset](#description-du-dataset)
- [Avant de Commencer le Travail](#avant-de-commencer-le-travail)
    - [Créer un environnement virtuel python](#creer-un-environnement-virtuel-pythoncreer-)
    - [Lancer L'environnement Python](#lancer-lenvironnement-python)
    - [Installer jupyter et jupyter kernel](#installer-jupyter-et-jupyter-kernel)
    - [Créer et installer un nouveau kernel jupyter](#creer-et-installer-un-nouveau-kernel-jupyter)
- [Lancer le code source](#lancer-le-code-source)
- [Resultats obtenus](#resultats-obtenus)
- [Références](#references)
## Description du dataset
INFORMATIONS SUR LES COLONNES LE DATASET
* stag : Expérience (temps)
* event : Rotation du personnel (démission de l'employé)
* gender : Sexe de l'employé (femme (f) ou homme (m))
* age : Âge de l'employé (année)
* industry : Industrie de l'employé
* profession : Profession de l'employé

## Avant de Commencer le Travail
Installet matlab et ajouter matlab au PATH
### Conversion du fichier csv en .mat
```matlab -nodesktop -r  ConverCSVTOMAP,quit```
ou
```matlab -nodesktop -r  Conver datasets/Saison20232024.csv Saison20232024 ,quit```
### Installer jupyter et jupyter kernel
```
conda install jupyter
conda install ipython
conda install ipykernel
pip install bash_kernel
```
### Créer et installer un nouveau kernel jupyter 
```
ipython kernel install --user --name=tp_map6014
python -m ipykernel install --user --name=tp_map6014
python -m bash_kernel.install
```
### Setup de l'environement suivant votre system
font et librairies
#### Ubuntu  
`sh setup/ubuntu/setup.sh`
#### Mac os
`bash setup/macos/setup.sh` 
## Lancer le code source (jupyter)
`Make open`
## Resultats obtenus
![courbe roc](./outputs/images/png/roc_curve.png)
## Avant le build
Installation des dependances du code python
`Make dep`
## Build
Vous pouvez build le model, pdf...
`Make build`
## Références
* [venv installations and kernel creation](https://medium.com/@WamiqRaza/how-to-create-virtual-environment-jupyter-kernel-python-6836b50f4bf4)
* [Theme](./theme/AdobeColor-My%20Color%20Theme-3.jpeg)
* [Latex](https://guides.nyu.edu/LaTeX/installation)
* [Color Generator](https://color.adobe.com)