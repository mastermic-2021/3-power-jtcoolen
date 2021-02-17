# Power

Un message a été encodé lettre à lettre dans une matrice `m`
d'éléments de Z/27Z, avec le codage devenu habituel (espace=0,a=1...).
Le chiffré correspond à la matrice `m^65537`. Son contenu se trouve dans
le fichier `input.txt`.


## Format

Vous devez écrire un programme `main.gp` en Pari/GP dont l'exécution via
```
gp -fq < main.gp
```
affiche (uniquement) la solution cherchée.

taper `make check` permet de vérifier que votre solution est bonne.

Veillez à mettre des commentaires sur votre démarche et sa validité
dans le fichier ``main.gp``.

De manière alternative, vous pouvez écrire un programme `main.c` qui
fasse la même chose.

## Validation

Il reste à faire un commit et à envoyer votre solution pour l'enregistrer
```
git add main.gp
git commit -m 'ma solution'
git push
```

