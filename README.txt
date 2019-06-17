

Ce programme utilise Pythia8 pour simuler des événements. Il faut déterminer les paramètres des événements grâce aux fichiers "file.cmnd" dans le dossier FILES;
Par défaut, 50 000 evenements sont simulés;
Plusieurs Macro root sont fournis, leur stucture global est fonctionnelle, il faut cependant les modifier pour qu'ils correspondent à ce que l'utilisateur souhaite représenter.



 3 cas sont possibles pour ce programme :

 							Cas 1 :Etude d’une seule simulation.

Dans ce cas, il n’y aura qu’une seule simulation effectuée qui piochera ses paramétres dans le
fichier "file.cmnd".
L’intégrale des données des particules est stockée dans le TTree "event" du fichier "Simu.root".
Les histogrammes des énergies sont quant à eux stockés dans Results.root; C’est le cas le plus
rapide pour avoir une vision du comportement à une certaine énergie.
De plus, ce cas utilise ses propres fichiers "Simu.root" et "file.cmnd".

							Cas 2: Etude de plusieurs simulations.

Il s’agit ici d’effectuer plusieurs simulation à différentes énergies dans le centre de masse afin
d’étudier la résolution totale en fonction de l’énergie.
Dans ce cas, il y aura autant de simulations que de fichiers "fileEcdm.cmnd" dans la liste au
début du main (par défaut 14, modifiable).
Le fichier "SimuEcdm.root" dans le dossier "Simulations" contiendra le TTree de la simulation
effectuée à l’énergie Ecdm.
Tout les histogrammes des énergies de chaque simulation sont dans le fichier "Results.root",
numérotés dans l’ordre. On peut donc dans ce cas récupérer des données à chaque énergie et
ainsi les tracer en fonction de celle-ci avec l'aide de MacroED2.C.

					Cas 3: Etude de la variation de la résolution totale en fonction d’un paramétre.

Ici, il va s’agir pour chaque simulation effectuée de regarder les effets de la variation d’un
paramétre de résolution sur la résolution totale.
Pour cela, l’utilisateur dispose d’une variable ResolVar qu’il pourra placer où bon lui semble
dans les formules de résolution et choisir les variation dans la boucle de la fonction "PlotEn-
ergy()".
Par défaut, ce test est effectué pour 3 énergies différentes (200 GeV, 350 Gev et 500 Gev).
On peut tracer ensuite en récupérant les 3 fichiers output les graphiques avec MacroED4.C.


Il existe un document qui explique exactement tous les détails du code et comment il fonctionne.



# Stage-2019-M1-ILC
