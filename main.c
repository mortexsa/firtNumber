#include <gmp.h>
#include <stdio.h>
#include <time.h>

void pgcd (mpz_t resultat, mpz_t a, mpz_t b) {
	mpz_t r;
	mpz_init(r);
	while (mpz_cmp_ui(b,0) != 0) {
		mpz_mod (r, a, b);
		mpz_set(a,b);
		mpz_set(b,r);
	}
	mpz_clear(r);
	mpz_set(resultat, a);
}

void exponentiationModulaire (mpz_t resultat, mpz_t a, mpz_t k, mpz_t n) {
	mpz_t i,tmp,tmpa;
	mpz_inits(i,tmp,tmpa,NULL);
	mpz_set(tmpa,a);
	for (mpz_set_ui(i,1); mpz_cmp_ui(k,0)>0; mpz_div_ui(k,k,2)) {
		mpz_mod_ui(tmp,k,2);
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

void jacobiSymbol (mpz_t resultat, mpz_t a, mpz_t b) {
	mpz_t tmp,tmp2,i;
	mpz_inits(tmp,tmp2,i,NULL);
	mpz_mod_ui(tmp,b,2);
	if ( mpz_cmp_ui(b,0) <= 0 || mpz_cmp_ui(tmp,0) == 0 )
		mpz_set_ui(resultat, 0);
	mpz_set_ui(i,1);
	if( mpz_cmp_ui(a,0) < 0 ) {
		mpz_neg( a, a);
		mpz_mod_ui(tmp,b,4);
		if (( mpz_cmp_ui(tmp,3) == 0)) 
			mpz_neg(i,i);
	}
	while (mpz_cmp_ui(a,0) != 0) {
		mpz_mod_ui(tmp,a,2);
		while (mpz_cmp_ui(tmp,0) == 0) {
			mpz_div_ui(a,a,2);
			mpz_mod_ui(tmp,b,8);
			if ((mpz_cmp_ui(tmp,3) == 0 || mpz_cmp_ui(tmp,5) == 0)) 
				mpz_neg(i,i);
			mpz_mod_ui(tmp,a,2);
		}
		mpz_set(tmp,a);
		mpz_set(a,b);
		mpz_set(b,tmp);
		mpz_mod_ui(tmp,a,4);
		mpz_mod_ui(tmp2,b,4);
		if (mpz_cmp_ui(tmp,3) == 0 && mpz_cmp_ui(tmp2,3) == 0)
			mpz_neg(i,i);
		mpz_mod(a,a,b);
	}
	if (mpz_cmp_ui(b,1) == 0){ mpz_set(resultat,i); }
	else mpz_set_ui(resultat,0);
	mpz_clears(tmp,tmp2,i,NULL);
}

int solovayStrassen(mpz_t aTraiter, mpz_t iterations) {
	
	if (mpz_cmp_ui(aTraiter,2) < 0)
		return 0;
	
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mod_ui(tmp,aTraiter,2);
	if (mpz_cmp_ui(tmp,0) == 0){
		mpz_clear(tmp);
		return 0;
	}

	mpz_t i,randomNumber,resultatJ,exposant,resultatM,randomTmp,aTraiterTmp;
	mpz_inits(i,randomNumber,resultatJ,exposant,resultatM,randomTmp,aTraiterTmp,NULL);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));

	mpz_sub_ui(exposant,aTraiter,1);
	mpz_div_ui(exposant,exposant,2);

	mpz_sub_ui(tmp, aTraiter, 3);

	for (mpz_set_ui(i,0); mpz_cmp(i,iterations) < 0; mpz_add_ui(i, i, 1)) {
		mpz_urandomm(randomNumber,state,tmp);
		mpz_add_ui(randomNumber, randomNumber, 2);
		mpz_set(aTraiterTmp,aTraiter);
		mpz_set(randomTmp,randomNumber);
		
		jacobiSymbol(resultatJ, randomTmp, aTraiterTmp);
		exponentiationModulaire(resultatM, randomNumber, exposant, aTraiterTmp);
		if(mpz_cmp_ui(resultatJ,0) == 0 && mpz_cmp(resultatJ,resultatM) != 0){

			mpz_clears(i,randomNumber,tmp,resultatJ,resultatM,randomTmp,aTraiterTmp,NULL);
			return 0;
		}
	}
	mpz_clears(i,randomNumber,tmp,resultatJ,resultatM,randomTmp,aTraiterTmp,NULL);
	return 1;
}

int main(){
	mpz_t r,a,b,res;
	mpz_inits(r,a,b,res,NULL);
	mpz_set_ui(a,33);
	mpz_set_ui(b, 100);
	// jacobiSymbol(r,a,b);
	// pgcd(r,a,b);
	// gmp_printf("%Zd \n",r);
	// exponentiationModulaire(res,r,a,b);
	// exponentiationModulaire(r,a,b);
	printf("%d \n",solovayStrassen(a,b));
	mpz_clears(r,a,b,res,NULL);


	return 0;
}
