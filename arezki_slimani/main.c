#include <gmp.h>
#include <stdio.h>
#include <time.h>

/**
* Cette Fonction permet de calculer un nombre 'a' a la puissance 'b'
* @param resultat : permet de renvoyer le resultat.
* @param nombre 
* @param exposant : l'exposant du nombre.
*/
void puissanceMoinUn(mpz_t resultat, mpz_t nombre, mpz_t exposant) {
	mpz_t i;
	mpz_init(i);
	mpz_set(resultat, nombre);
	for (mpz_set_ui(i,1); mpz_cmp(i,exposant)<0; mpz_add_ui(i,i,1)) {
		//On le multiplie par lui meme 'exposant fois'.
		mpz_mul(resultat,resultat,nombre);
	}
	mpz_sub_ui(resultat,resultat,1);
	mpz_clear(i);
}

/**
* Cette fonction permet de Calculer l'exponentiation Modulaire [a^k mod n].
* @param resultat : permet de renvoyer le resultat.
* @param a : le nombre qu'on mettra a la puissance k.
* @param k : exposant du nombre a
* @param n : on met le tout modulo n.
*/

void exponentiationModulaire (mpz_t resultat, mpz_t a, mpz_t k, mpz_t n) {
	mpz_t i,tmp,tmpa;
	mpz_inits(i,tmp,tmpa,NULL);
	mpz_set(tmpa,a);
	//parcours de boucle k = k/2
	for (mpz_set_ui(i,1); mpz_cmp_ui(k,0)>0; mpz_div_ui(k,k,2)) {
		mpz_mod_ui(tmp,k,2);
		// Si l'exposant k est different de 0
		if(mpz_cmp_ui(tmp,0) != 0){
			mpz_mul(i,i,tmpa);
			mpz_mod(i,i,n);
		}
		mpz_mul(tmpa,tmpa,tmpa);
		mpz_mod(tmpa,tmpa,n);
	}
	mpz_set(resultat,i);
	mpz_clears(i,tmp,tmpa,NULL);
}

/**
* Cette fonction permet de calculer le Symbole de jacobi (a/p) et de determiner si
* si p divise a ou pas ainsi que si a est un résidu quadratique modulo p ou non
* @param resultat : On renvoi le resultat.
* @param a : a est un résidu quadratique ou non de b
* @param b 
*/
void jacobiSymbol (mpz_t resultat, mpz_t a, mpz_t b) {
	mpz_t tmp,tmp2,i;
	mpz_inits(tmp,tmp2,i,NULL); // on initialise la liste de variable
	mpz_mod_ui(tmp,b,2); // tmp = b % 2
	//Si b % 2 est egale a 0 ou b = 0
	if ( mpz_cmp_ui(b,0) <= 0 || mpz_cmp_ui(tmp,0) == 0 )
		mpz_set_ui(resultat, 0);
	mpz_set_ui(i,1);
	if( mpz_cmp_ui(a,0) < 0 ) { //Si a < 0
		mpz_neg( a, a);
		mpz_mod_ui(tmp,b,4);
		if (( mpz_cmp_ui(tmp,3) == 0)) //si b % 4 = 3
			mpz_neg(i,i); // i =- i
	}
	while (mpz_cmp_ui(a,0) != 0) {
		mpz_mod_ui(tmp,a,2);
		while (mpz_cmp_ui(tmp,0) == 0) { //Tant que a % 2 = 0
			mpz_div_ui(a,a,2);
			mpz_mod_ui(tmp,b,8);
			if ((mpz_cmp_ui(tmp,3) == 0 || mpz_cmp_ui(tmp,5) == 0)) //Si b % 8 = 3 ou b % 8 = 5
				mpz_neg(i,i);
			mpz_mod_ui(tmp,a,2);
		}
		mpz_set(tmp,a);
		mpz_set(a,b);
		mpz_set(b,tmp);
		mpz_mod_ui(tmp,a,4);
		mpz_mod_ui(tmp2,b,4);
		if (mpz_cmp_ui(tmp,3) == 0 && mpz_cmp_ui(tmp2,3) == 0) //Si a % 4 = 3 et b % 4 = 3
			mpz_neg(i,i);
		mpz_mod(a,a,b);
	}
	if (mpz_cmp_ui(b,1) == 0){ mpz_set(resultat,i); } //Si b = 1
	else mpz_set_ui(resultat,0);
	mpz_clears(tmp,tmp2,i,NULL); // on libere la mémoire 
}

/**
* Cette fonction permet de calculer si un nombre est premier ou composé.
* @param aTraiter : le nombre que l'on désire traiter.
* @param iterations : le nombre d'iterations que l'on désire faire lors du teste.
* @return on return 1 s'il est premier ou bien 0 s'il est composé
*/
int solovayStrassen(mpz_t aTraiter, mpz_t iterations) {
	
	if (mpz_cmp_ui(aTraiter,2) < 0)
		return 0;
	
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mod_ui(tmp,aTraiter,2);
	if (mpz_cmp_ui(tmp,0) == 0){ // SI aTraiter % 2 est egale a 0
		mpz_clear(tmp);
		return 0;
	}

	mpz_t i,randomNumber,resultatJ,exposant,resultatM,randomTmp,aTraiterTmp;
	mpz_inits(i,randomNumber,resultatJ,exposant,resultatM,randomTmp,aTraiterTmp,NULL);
	//on initialise tout les parametre pour avoir des nombres aleatoire
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));

	mpz_sub_ui(exposant,aTraiter,1); //On fait exposant-1
	mpz_div_ui(exposant,exposant,2); // puis on fait (exposant-1)/2

	mpz_sub_ui(tmp, aTraiter, 2); // on soustrait 2 pour avoir l'ensemble de définition compris entre 2 et n-1 pour les nombre aléatoire

	for (mpz_set_ui(i,0); mpz_cmp(i,iterations) < 0; mpz_add_ui(i, i, 1)) {
		//on creer notre nombre aleatoire
		mpz_urandomm(randomNumber,state,tmp);
		mpz_add_ui(randomNumber, randomNumber, 2); // on reajoute le 2 qu'on a soustrait precedemment
		//On utilise des variable de copie
		mpz_set(aTraiterTmp,aTraiter);
		mpz_set(randomTmp,randomNumber);
		jacobiSymbol(resultatJ, randomTmp, aTraiterTmp);
		exponentiationModulaire(resultatM, randomNumber, exposant, aTraiterTmp);
		// Si jacobie donne 0 et que jacobi est different de l'exponentiation modulaire, alors on renvoi 0
		if(mpz_cmp_ui(resultatJ,0) == 0 && mpz_cmp(resultatJ,resultatM) != 0){

			mpz_clears(i,randomNumber,tmp,resultatJ,resultatM,randomTmp,aTraiterTmp,NULL);
			return 0;
		}
	}
	mpz_clears(i,randomNumber,tmp,resultatJ,resultatM,randomTmp,aTraiterTmp,NULL);
	return 1;
}



int main(){
	//Ceci concerne l'affichage sur le terminal.
	mpz_t recup,iterations,exposant,deux,aAnalyser;
	mpz_inits(recup,iterations,exposant,deux,aAnalyser,NULL);
	gmp_printf("////////////// BIENVENU AU TESTEUR DE NOMBRES PREMIERS //////////////\n");
	gmp_printf("Tapez : \n  1 Pour acceder a la liste de teste deja présent.\n  2 Pour entrer le nombre de votre choix.\n");
	do {
		mpz_set_ui(recup,0);
		gmp_printf("Je choisi : ");
		gmp_scanf("%Zd", recup);
	}while (mpz_cmp_ui(recup,1) != 0 && mpz_cmp_ui(recup,2) != 0);
	if(mpz_cmp_ui(recup,1) == 0){
		gmp_printf("Voici la liste des 10 nombres disponibles,Tapez le chiffre qui lui correspend pour le choisir :\n");
		gmp_printf("  1-[2^521-1]\n  2-[2^1279-1]\n  3-[2^2281-1]\n  4-[2^3217-1]\n  5-[2^4423-1]\n  6-[2^9941-1]\n  7-[2^11213-1]\n  8-[2^19937-1]\n  9-[2^21701-1]\n  10-[2^23209-1]\n");
		
		do {
			mpz_set_ui(recup,0);
			gmp_printf("Entrez le numero: ");
			gmp_scanf("%Zd",recup);
		}while((mpz_cmp_ui(recup,11)) >= 0 || (mpz_cmp_ui(recup,1)) < 0);
		gmp_printf("Entrez le nombre d'iterations : ");
		gmp_scanf("%Zd",iterations);
		
			if(mpz_cmp_ui(recup,1)==0)
				mpz_set_ui(exposant,521);
			else if(mpz_cmp_ui(recup,2)==0)
				mpz_set_ui(exposant,1279);
			else if(mpz_cmp_ui(recup,3)==0)
				mpz_set_ui(exposant,2281);
			else if(mpz_cmp_ui(recup,4)==0)
				mpz_set_ui(exposant,3217);
			else if(mpz_cmp_ui(recup,5)==0)
				mpz_set_ui(exposant,4423);
			else if(mpz_cmp_ui(recup,6)==0)
				mpz_set_ui(exposant,9941);
			else if(mpz_cmp_ui(recup,7)==0)
				mpz_set_ui(exposant,11213);
			else if(mpz_cmp_ui(recup,8)==0)
				mpz_set_ui(exposant,19937);
			else if(mpz_cmp_ui(recup,9)==0)
				mpz_set_ui(exposant,21701);
			else if(mpz_cmp_ui(recup,10)==0)
				mpz_set_ui(exposant,23209);
		
		mpz_set_ui(deux,2);
		puissanceMoinUn(aAnalyser,deux,exposant);
	}
	else {
		gmp_printf("Entrez le nombre de votre choix : \n");
		gmp_printf("Le nombre est :");
		gmp_scanf("%Zd",aAnalyser);
		gmp_printf("Nombre d'iterations :");
		gmp_scanf("%Zd",iterations);
	}
	int leResultat = 0;
	leResultat = solovayStrassen(aAnalyser,iterations);
	if(leResultat == 1){
		gmp_printf("Ce nombre est probablement Premier.\n");
	}
	else {
		gmp_printf("Ce nombre est Composé.\n");
	}

	mpz_clears(recup,iterations,exposant,deux,aAnalyser,NULL);
	return 0;
}
