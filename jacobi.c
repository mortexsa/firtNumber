#include <gmp.h>
#include <stdio.h>

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
			mpz_divexact_ui(a,a,2);
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

int main(){
	mpz_t r,a,b;
	mpz_inits(r,a,b,NULL);
	mpz_set_ui(a, 10);
	mpz_set_ui(b, 5);
	jacobiSymbol(r,a,b);
	gmp_printf("%Zd \n",r);
	mpz_clears(r,a,b,NULL);


	return 0;
}