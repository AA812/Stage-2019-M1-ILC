

Ce programme utilise Pythia8 pour simuler des événements. Il faut déterminer les paramètres des événements grâce aux fichiers "file.cmnd" dans le dossier FILES;
Par défaut, 50 000 evenements sont simulés;
Plusieurs Macro root sont fournis, leur stucture global est fonctionnelle, il faut cependant les modifier pour qu'ils correspondent à ce que l'utilisateur souhaite représenter.





 3 cas sont possibles pour ce programme :

 	1. Etude d'une seule simulation :

   		Dans ce cas, il n'y aura qu'une seule simulation effectué qui piochera ses paramétres dans le fichier "file.cmnd";
    		L'intégrale des données des particules est stocké dans le TTree "event" du fichier Simulation.root;
   		Les histogrammes des énergies sont quant à eux stockés dans Histos.root;


	 2. Etude de plusieurs simulations : 
 
    		Il s'agit ici d'effectuer plusieurs simulation à différentes énergies dans le centre de masse afin d'étudier la résolution totale en fonction 			de l'énergie;
    		Dans ce cas, il y aura autant de simulation que de fichier "fileEcdm.cmnd" dans la liste au début du main (par défaut 14, modifiable);
                Le fichier Simulation.root ne contiendra que le TTree de la dernière simulation effectuée;
                Cependant tout les histogrammes de chaque simulation sont dans le fichier Histos.root, numérotés dans l'ordre.
    		Le fichier "Output.txt" sert à stocker les points qui serviront à tracer le graph ResTot=f(E);
    		Pour tracer celui-ci, une fois le programme terminé, éxécuter "MacroED2.C" dans root;


 	3. Etude de la variation de la résolution totale en fonction d'un paramétre :

   	 	Ici, il va s'agir pour chaque simulation effectuée de regarder les effets de la variation d'un paramétre de résolution sur la résolution 			totale;
    		Pour cela, l'utilisateur dispose d'une variable ResolVar qu'il pourra placer où bon lui semble dans les formules de résolution et choisir les 			variation dans la boucle de la fonction "PlotEnergy();
   		Par défaut, ce test est effectué pour 3 énergies différentes, les valeurs étant stockées dans 3 fichiers différents : "Output1,2,3.txt";
    		Pour tracer ResTot=f(ResolVar) éxécuter MacroED4.C dans root;
# Stage-2019-M1-ILC
