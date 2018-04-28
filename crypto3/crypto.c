#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

//PGCD
void pgcd(mpz_t result,mpz_t a, mpz_t b) // result=resultat, a=entier, b=entier

{
		mpz_mod(result,a,b);//result= a%b
		
		if(mpz_cmp_ui(result,0)==0) // Si result == 0
		{
			mpz_set(result,b); //return b
		}
		else
		{
			pgcd(b,b,result); // Recursif sur pgcd(b,result)
		}
}

//Legendre
void legendre(mpz_t result,mpz_t a,mpz_t p)
{
	mpz_t reste,p1,a1;
	mpz_inits(reste,p1,a1,NULL);
	
	 if(mpz_cmp_ui(p,2)<0) //Si p<2
	 {
		 gmp_printf("Erreur, p ne doit pas être inférieur a 2");//Erreur
	 }
	 
	 if(mpz_cmp_ui(a,0)==0 || mpz_cmp_ui(a,0)==1)//Si a=0 ou a=1
	 {
		 mpz_set(result,a);// result = a
	 }
	 
		mpz_mod_ui(reste,a,2);// reste= a % 2
		//gmp_printf("modulo a =%Zd \n",reste);
	
	 if(mpz_cmp_ui(reste,0)==0) // if a % 2 =0
	 {
		 mpz_div_ui(a,a,2);// a=a/2
		 legendre(result,a,p);//result=legendre(a/2,p)
		  //~ gmp_printf("je suis dans le bendo");
		 mpz_mul(p1,p,p);//p1=p*p
		 mpz_sub_ui(p1,p1,1);//p1=p1-1
		 
		 if(mpz_cmp_ui(p1,8)==0) 
		 {
			 gmp_printf("je suis dans le bendo");
			 mpz_neg(result,result);
		 }
	 } 
	 
	 else
	 {
		 mpz_mod(reste,p,a);//reste=p % a 
		 gmp_printf("reste=%Zd",reste);
		 //~ legendre(result,reste,a);//result=legendre(p%a,a)
		 
		 //~ mpz_sub_ui(a1,a,1);//a1=a-1
		 //~ mpz_sub_ui(p1,p,1);//p1=p-1
		 //~ mpz_mul(a1,a1,p1);//a1=a1*p1
		 
		 //~ if(mpz_cmp_ui(a1,4)==0)
		 //~ {
			 //~ mpz_neg(result,result);
		 //~ }
	 }
	
}

//Exponentiation Modulaire
void deux_pow_exp(mpz_t result, mpz_t exp){ // Calcul 2^exp et stock le résultat dans result. 
 	 
 	mpz_t i; 
 	mpz_init(i); 
	mpz_set_ui(result, 1); // Result = 1  
 
 
 	for (mpz_set_ui(i,1); mpz_cmp(i,exp) <= 0; mpz_add_ui(i, i, 1))//for (i= 1; i <= exp; i++) 
 	{ 
 		mpz_mul_ui(result,result,2); // result = result*2 
 	}
 	 
	mpz_clear(i); 
} 

//Square And Multiply
void square_multiply(mpz_t result, mpz_t a, mpz_t n, mpz_t h){ // result=resultat de la fonction, a=entier , n=modulo, h=exposant  

	mpz_t reste;	// "reste" va contenir le reste de la division euclidienne de exp par 2 à chaque itération.	 
 	mpz_init(reste); 
 	mpz_set(result,a);	//result = a 
 
 	int i = 0; 
	int* tab = malloc(sizeof(int)*1024); // allocation mémoire d'un tableau de taille 1024 
	// Convertion de h en binaire 
	while( mpz_cmp_ui(h,0) > 0 ){ 
 		mpz_fdiv_qr_ui(h, reste, h, 2); // reste = reste de la div et h = h/2 pour la prochaine itération 
 		if(mpz_cmp_ui(reste,0) == 0){ 
			tab[i] = 0; 
		}else{ 
			tab[i] = 1; 
 		} 
		i++; 
	} 
		i--; // on ne prend pas en compte t  
 
 	while (i > 0){ 	//mpz_cmp_ui(exp, 1) > 0){ 
 		mpz_mul(result, result, result); // result = result*result 
		mpz_mod(result, result, n); // result = result*result modulo n 
		if ( tab[--i] == 1 ){ // Bit à 1 
			mpz_mul(result, result, a); // result = result*a 
			mpz_mod(result, result, n); // result = result*a modulo n 
		} 
 	} 
	 
 	// Libération mémoire 
	free(tab); 
 	tab = NULL; 
	mpz_clear(reste); 
} 
void change(mpz_t a,mpz_t b){
	mpz_t tmp;
	mpz_init(tmp);
	mpz_set(tmp,a);
	mpz_set(a,b);
	mpz_set(b,tmp);
}

int jacobiversion2(mpz_t a,mpz_t b)
{
   mpz_t tmp;
   mpz_t tmp2;
   mpz_inits(tmp,tmp2,NULL);
   mpz_mod_ui(tmp,b,2);
   
    if(mpz_cmp_ui(b,0)<=0 || mpz_cmp_ui(tmp,0)==0) return 0;
   mpz_t j;
   mpz_init(j);
   mpz_set_ui(j,1);
   
   if (mpz_cmp_ui(a,0)==0) {
	  mpz_neg(a,a);
	  mpz_mod_ui(tmp,b,4);
      if ( mpz_cmp_ui(tmp,3)==0) mpz_neg(j,j);
      }
  
  while(mpz_cmp_ui(a,0)!=0) {
      while(mpz_cmp_ui(tmp,0)==0) {
         /* Process factors of 2: Jacobi(2,b)=-1 if b=3,5 (mod 8) */
         mpz_div_ui(a,a,2);
         mpz_mod_ui(tmp,b,8);
         if( mpz_cmp_ui(tmp,3)==0 || mpz_cmp_ui(tmp,5)==0) mpz_neg(j,j);
         mpz_mod_ui(tmp2,a,2);//tmp2=a%2
         }
      /* Quadratic reciprocity: Jacobi(a,b)=-Jacobi(b,a) if a=3,b=3 (mod 4) */
      change(a,b);
      mpz_mod_ui(tmp,a,4);//tmp=a%4
      mpz_mod_ui(tmp2,b,4);//tmp=a%4

      if(mpz_cmp_ui(tmp,3)==0 && mpz_cmp_ui(tmp2,3)==0)  mpz_neg(j,j);
      mpz_mod(a,a,b);
							}
   if (mpz_cmp_ui(b,1)==0) { gmp_printf("Jacobi=%Zd",j); return -1;}
gmp_printf("Jacobi=0");
return 0;
}

int jacobi(mpz_t a,mpz_t n)
{
	if(mpz_cmp_ui(a,0)==0)
	{
		return 0;
	}
	
	mpz_t result,temp,reste,reste1,reste2,reste3,reste4,reste5;
	mpz_inits(result,temp,reste,reste1,reste2,reste3,reste4,reste5,NULL);
	mpz_set_ui(result,1);
	
	if(mpz_cmp_ui(a,0)<0)
	{
		mpz_neg(a,a);
		
		mpz_mod_ui(reste,n,4);
		
		if(mpz_cmp_ui(reste,0)== 3)
		{
			mpz_neg(result,result);
			return -1;
		}
	}
	
	if(mpz_cmp_ui(a,1)==0)//Propriété 4
	{
		mpz_set(result,result);
		return 1;
	}

	
	while(mpz_cmp_ui(a,1)==0)
	{
		if(mpz_cmp_ui(a,0)<0)
		{
			mpz_neg(a,a);
			
			mpz_mod_ui(reste,n,4);
			
			if(mpz_cmp_ui(reste,0)== 3)
			{
				mpz_neg(a,a);
			}
		
		}
		
		mpz_mod_ui(reste,a,2);
		
		while(mpz_cmp_ui(reste,0)==0)
		{
			mpz_div_ui(a,a,2);
			mpz_mod_ui(reste1,n,8);
			mpz_mod_ui(reste2,n,5);

			if(mpz_cmp_ui(reste1,3)==0 || mpz_cmp_ui(reste2,5)==0)
			{
				mpz_neg(result,result);
			}
			
			change(a,n);
			
			mpz_mod_ui(reste3,a,4);
			mpz_mod_ui(reste4,n,4);
			
			if(mpz_cmp_ui(reste3,3)==0 || mpz_cmp_ui(reste4,3)==0)
			{
				mpz_neg(result,result);
			}
			
				mpz_mod(reste5,a,n);
				mpz_set(a,reste5);
				
				mpz_div_ui(n,n,2);//n=n/2
				
			if(mpz_cmp(a,n)>0)
			{
				mpz_sub(a,a,n);
			}
		}

	}
	
		if(mpz_cmp_ui(n,1)==0)
		{
			mpz_neg(result,result);
		}
}

int solovay_strassen(mpz_t n,mpz_t iteration)
{
	mpz_t reste,i,a,iteration_3,val_sm,result_j,result_sm,n1,n2;
	mpz_inits(reste,i,a,iteration_3,val_sm,result_j,result_sm,n1,n2,NULL);
	
	if(mpz_cmp_ui(n,2)<0)
	{
		return -1;
	}
	
	mpz_mod_ui(reste,n,2);
	
	if((mpz_cmp_ui(n,2)!=0) && (mpz_cmp_ui(reste,0)==0))
	{
		return -1;
	}
	
	gmp_randstate_t etat; 
 	gmp_randinit_mt(etat); 
 	gmp_randseed_ui(etat,time(NULL));
 	
 	
 	mpz_sub_ui(iteration_3,iteration,3);//iteration_1=iteration-3

	
	for(mpz_set_ui(i,0); mpz_cmp(i,iteration) <= 0; mpz_add_ui(i, i, 1))//for (i= 0; i < iteration; i++)
	{
		mpz_urandomm(a,etat,iteration_3); // 0 <= a <= iteration-3
		mpz_add_ui(a,a,2);//a+2
		
		int j=jacobi(a,n);//jacobi(a,n);
		mpz_set_ui(result_j,j);
		
		mpz_sub_ui(n1,n,1);//n1=n-1
		mpz_div_ui(n2,n1,2);//n2=n1/2
		
		square_multiply(val_sm,a,n,n2);
		mpz_set(result_sm,val_sm);
		
		
		if(mpz_cmp_ui(result_j,0)==0 || mpz_cmp(result_j,result_sm)!=0)
		{
			return 0;
		}
		
		else
		{
			return 1;
		}
	}
	
	return 1;
	
}
   

int main(int argc, char **argv)
{
	//Test PGCD
	mpz_t a;
	mpz_t b;
	mpz_t r;
	
	mpz_inits (a,b,r,NULL);

	mpz_set_ui(a,1234);
	mpz_set_ui(b,357);
	
	gmp_printf("PGCD(%Zd,%Zd)",a,b);
	
	pgcd(r,a,b);
	
	gmp_printf("=%Zd \n",r);
	
	mpz_clears (a,b,r,NULL);
	
	//Fin Test PGCD OK
	
	//Test 2^exp
	
	mpz_t exp;
	mpz_t r2;
	
	mpz_inits (exp,r2,NULL);
	
	mpz_set_ui(exp,521);
	
	deux_pow_exp(r2,exp);
	
	mpz_sub_ui(r2,r2,1);//r2=2^exp -1
	
	gmp_printf("2^%Zd-1",exp);
	
	gmp_printf("=%Zd\n",r2);
	
	
	
	//Fin Test 2^exp OK
	
	//Test Legendre
	
	mpz_t r3;
	mpz_t c;
	mpz_t p;
	
	mpz_inits (r3,c,p,NULL);

	mpz_set_ui(c,10);
	mpz_set_ui(p,3);
	
	legendre(r3,c,p);
	
	gmp_printf("Legendre=%Zd\n",r3);
	
	mpz_clears(r3,c,p,NULL);

	
	//Fin Test Legendre
	
	//Test jacobi
	mpz_t A;
	mpz_t B;
	mpz_t r4,j,h,g;
	
	mpz_inits (A,B,r4,j,h,g,NULL);
	
	mpz_set_ui(A,8721);
	mpz_set_ui(B,4235);
	mpz_set_ui(j,3);
	mpz_set_ui(h,143);
	mpz_set_ui(g,103);

//~ int okok=mpz_jacobi(A,B);
//~ printf("jacobi=%d",okok);
	
		int okok=jacobi(A,B);
		printf("Jacobi=%d",okok);
		
	
		
	//~ square_multiply(r4,j,h,g);
	//~ gmp_printf("/n%Zd",r4);
	
	//~ mpz_clears(A,B,j,h,g);
     
     //Test Solovay-Strassen
     
     int final;
     mpz_t n,iteration,ok,k;
     mpz_inits(n,ok,k,NULL);
     
     mpz_set_ui(ok,7);
     gmp_printf("\n%Zd=",ok);
     mpz_set_ui(k,4000);
     final=solovay_strassen(r2,k);//r2=(2^521)-1
     
     if(final==1)
     {
		 gmp_printf("Ce nombre est Premier");
	 }
	 
	 if(final==0)
	 {
		 gmp_printf("Ce nombre est Composé");
	 }
	 
	 if(final==-1)
	 {
		 gmp_printf("\nErreur");
	 }
     
     mpz_clears (exp,r2,k,NULL);
     
     //Fun Solovay Strassen

	
	return 0;
}

