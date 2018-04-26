#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

//PGCD
void pgcd(mpz_t result,mpz_t a,mpz_t b)//Calcul Du PGCD 
{
	mpz_t a2,b2;//Declaration a2,b2 (Variable temporaire)
	mpz_inits(a2,b2,NULL);//Initialisation a2,b2
	
	mpz_set(a2,a);//a2=a
	mpz_set(b2,b);//b2=b

	while(mpz_cmp_ui(b2,0)==1)
	{
		mpz_mod(result,a2,b2);//result= a2%b2
		mpz_set(a2,b2);// a2=b2
		mpz_set(b2,result);//b2=result
	}

	mpz_set(result,a2);//result=a2
	
	mpz_clears(a2,b2,NULL);//Libération Memoire a2,b2
}

//Exposant
void deux_pow_exp(mpz_t result,mpz_t exp)// Calcul 2^exp et stock le résultat dans result.
{  
 	mpz_t i; //Declaration i
 	mpz_init(i);//Initialisation i
 	 
	mpz_set_ui(result,1); //result = 1  
 
 	for (mpz_set_ui(i,1); mpz_cmp(i,exp) <= 0; mpz_add_ui(i,i,1))//for (i= 1; i <= exp; i++) 
 	{ 
 		mpz_mul_ui(result,result,2); //result = result*2 
 	}
 	 
	mpz_clear(i); //Liberation Memoire i
}

//Symbole De Jacobi
void jacobi(mpz_t result,mpz_t a,mpz_t n)//Calcul Du Symbole De Jacobi
{
	mpz_t i,temp,temp2;//Declaration tmp,tmp_2,i
	mpz_inits(temp,temp2,i,NULL);//Initialisation tmp,tmp2,i
	
	mpz_set_ui(i,1);//i=1
	
	mpz_mod_ui(temp,n,2);//temp=n%2
	
	if (mpz_cmp_ui(n,0)<=0 || mpz_cmp_ui(temp,0)==0)//Si n<=0 OU temp ==0
	{
		mpz_set_ui(result,0);//result=0
	}
	
	if(mpz_cmp_ui(a,0)<0)//Si a<0
	{
		mpz_neg(a,a);//a=-a
		mpz_mod_ui(temp,n,4);//temp=n%4
		
		if((mpz_cmp_ui(temp,3)==0))//Si temp==3
		{ 
			mpz_neg(i,i);//i=-i
		}
	}
	
	while(mpz_cmp_ui(a,0)!=0)//Si a!=0
	{
		mpz_mod_ui(temp,a,2);//temp=a%2
		
		while(mpz_cmp_ui(temp,0)==0)//Tant que temp==0
		{
			mpz_divexact_ui(a,a,2);//a=a/2
			mpz_mod_ui(temp,n,8);//temp=n%8
			
			if((mpz_cmp_ui(temp,3) == 0 || mpz_cmp_ui(temp,5) == 0))//Si temp==3 OU temp==5
			{
				mpz_neg(i,i);//i=-i
			}
			
			mpz_mod_ui(temp,a,2);//temp=a%2
		}
		
		mpz_swap(a,n);//change a et b
		
		mpz_mod_ui(temp,a,4);//temp=a%4
		mpz_mod_ui(temp2,n,4);//temp2=n%4
		
		if (mpz_cmp_ui(temp,3) == 0 && mpz_cmp_ui(temp2,3) == 0)//Si temp==3 ET temp2==3
		{
			mpz_neg(i,i);//i=-i
		}
		
		mpz_mod(a,a,n);//a=a%n
	}
	
	if (mpz_cmp_ui(n,1) == 0)//Si n==1
	{ 
		mpz_set(result,i);//result=i
	}
	
	else 
	{
		mpz_set_ui(result,0);//result=0
	}
	
	mpz_clears(temp,temp2,i,NULL);//Liberation De Memoire
} 

//Square And Multiply
void square_multiply(mpz_t a,mpz_t n,mpz_t h)
{	
	mpz_t r,p;//Declaration r,p
	mpz_inits(r,p,NULL);//Initialisation r,p
	
	for(mpz_set_ui(p,1);mpz_cmp_ui(h,0)>0;mpz_div_ui(h,h,2))//for ( p=1;h>0,h=h/2)
	{
		mpz_mod_ui(r,h,2);//r=h%2
		
		if(mpz_cmp_ui(r,0)!=0)//Si r!=0
		{
			mpz_mul(p,p,a);//p=p*a
			mpz_mod(p,p,n);//p=p%n
		}
		
		mpz_mul(a,a,a);//a=a*a
		mpz_mod(a,a,n);//a=a%n
	}
	
	mpz_clears(r,p,NULL);//Liberation Memoire r,p
}

//Solovay-Strassen
int solovay_strassen(mpz_t n, mpz_t k)
{
	mpz_t a,a2,exp,i,r,n2,max,test;//Declaration des mpz
	mpz_inits(a,a2,i,r,exp,n2,max,test,NULL);//Initialisation des mpz
	
	mpz_sub_ui(exp,n,1);//exp=n-1
	mpz_div_ui(exp,n,2);//exp=(n-1)/2
	
	mpz_sub_ui(max,n,3);//n2=n-3
	mpz_mod_ui(test,n,2);
	
	if(mpz_cmp_ui(n,2)<0)
	{
		 return 0;
	}
	
	if( mpz_cmp_ui(test,0)==0) 
	{
		return 0;
	}

	gmp_randstate_t etat;//Declaration etat
	gmp_randinit_mt(etat);//Initialisation etat
	gmp_randseed_ui(etat,time(NULL));//etat=time(NULL)
      		
	for (mpz_set_ui(i,1);mpz_cmp(i,k)<=0;mpz_add_ui(i,i,1))//for (i=1;i<=k;i++)
	{
		mpz_urandomm(a,etat,max);// 0 <= a <= n-3
		mpz_add_ui(a,a,2);//2 <= a <= n-1	
			
		mpz_set(n2,n);//n2=n
		mpz_set(a2,a);//n2=n
		printf("....");
		jacobi(r,a2,n2);//jacobi
		square_multiply(a,n,exp);//square and multiply
		
		if(mpz_cmp_ui(r,0)==0 && mpz_cmp(a,r)!=0)
		{ 
				mpz_clears(n2,exp,max,a,a2,i,r,test,NULL);//Liberation Memoire
				return 0;
		}
	}
	
	mpz_clears(n2,exp,a,a2,i,max,r,test,NULL);//Liberation Memoire
	return 1; 	
}

//Affichage
void affichage()
{
	mpz_t exp,exp2,i,i2;
	mpz_inits(exp,exp2,i,i2,NULL);
	
	for (mpz_set_ui(i,1); mpz_cmp_ui(i,5) <= 0; mpz_add_ui(i, i, 1))//Affichage Choix Utilisateur	
	{
		if(mpz_cmp_ui(i,1)==0)
		{
			mpz_set_ui(exp,521);
			mpz_set_ui(exp2,9941);
			mpz_add_ui(i2,i,5);
		}
		
		if(mpz_cmp_ui(i,2)==0)
		{
			mpz_set_ui(exp,1279);
			mpz_set_ui(exp2,11213);
			mpz_add_ui(i2,i,5);
		}
		
		if(mpz_cmp_ui(i,3)==0)
		{
			mpz_set_ui(exp,2281);
			mpz_set_ui(exp2,19937);
			mpz_add_ui(i2,i,5);
		}
		
		if(mpz_cmp_ui(i,4)==0)
		{
			mpz_set_ui(exp,3217);
			mpz_set_ui(exp2,21701);
			mpz_add_ui(i2,i,5);
		}
		
		if(mpz_cmp_ui(i,5)==0)
		{
			mpz_set_ui(exp,4423);
			mpz_set_ui(exp2,23209);
			mpz_add_ui(i2,i,5);
		}
		
		gmp_printf("[%Zd]:(2^%Zd)-1     [%Zd]:(2^%Zd)-1\n",i,exp,i2,exp2);
	}
		gmp_printf("\nOU\n\n[11]:Nombre De Votre Choix\n");
		mpz_clears(exp,exp2,i,i2,NULL);
}

int main()
{
	int final=0;
	mpz_t nombre,res_jacobi,ja,r,jn,pgcda,pgcdb,exp,choix,exp2,i,i2,iteration,action;//Declaration nombre,res_jacobi,ja,jn,choix,exp,exp2,action,iteration,i,i2
	mpz_inits(nombre,res_jacobi,r,ja,jn,pgcda,pgcdb,choix,iteration,i,i2,exp,exp2,action,NULL);//Initialisation choix,res_jacobi,ja,jn,nombre,iteration,exp,exp2,action,i,i2
	
	gmp_printf("\n###################### Test PGCD ######################\n\n");
	
	mpz_set_ui(pgcda,1234);
	mpz_set_ui(pgcdb,357);
	
	pgcd(r,pgcda,pgcdb);
	
	gmp_printf("PGCD(%Zd,%Zd)=%Zd\n",pgcda,pgcdb,r);
	
	gmp_printf("\n###################### Test Jacobi ######################\n\n");
	
	mpz_set_ui(ja,8721);
	mpz_set_ui(jn,4235);

	gmp_printf("Jacobi(%Zd/%Zd)=",ja,jn);
	
	jacobi(res_jacobi,ja,jn);
	gmp_printf("%Zd\n",res_jacobi);
	
	mpz_set_ui(ja,541);
	mpz_set_ui(jn,2011);
	
	gmp_printf("Jacobi(%Zd/%Zd)=",ja,jn);
	
	jacobi(res_jacobi,ja,jn);
	gmp_printf("%Zd\n",res_jacobi);	
	
	mpz_set_ui(ja,8375);
	mpz_set_ui(jn,5371);
	
	gmp_printf("Jacobi(%Zd/%Zd)=",ja,jn);
	
	jacobi(res_jacobi,ja,jn);
	gmp_printf("%Zd\n",res_jacobi);
	
	mpz_set_ui(ja,0);
	mpz_set_ui(jn,5371);
	
	gmp_printf("Jacobi(%Zd/%Zd)=",ja,jn);
	
	jacobi(res_jacobi,ja,jn);
	gmp_printf("%Zd\n",res_jacobi);
	
	gmp_printf("\n###################### Test Programme ######################\n\n");
	
	mpz_set_ui(nombre,2);//nombre=2
	

	while(1)//Condition d'arret
	{		
		affichage();
		
		gmp_printf("\n**Choisir un numéro:");
		gmp_scanf("%Zd",action);//Demande à l'utilisateur le numéro
		
		if(mpz_cmp_ui(action,1000)==0)
		{
			 break;
		}
		
		if( mpz_cmp_ui(action,12)>=0 || mpz_cmp_ui(action,0)<=0) 
		{
			continue;
		}
		
		gmp_printf("**Choisir le nombre d'itérations:");
		gmp_scanf("%Zd",iteration);//Demande à l'utilisateur l'iteration
		
		if(mpz_cmp_ui(action,11)!=0)
		{
			if(mpz_cmp_ui(action,1)==0)//Choix N°1
			{
				mpz_set_ui(exp,521);
			}
					
			else if(mpz_cmp_ui(action,2)==0)//Choix N°2
			{
				mpz_set_ui(exp,1279);
			}
		
		else if(mpz_cmp_ui(action,3)==0)//Choix N°3
			{
				mpz_set_ui(exp,2281);
			}
		else if(mpz_cmp_ui(action,4)==0)//Choix N°4
			{
				mpz_set_ui(exp,3217);
			}
		else if(mpz_cmp_ui(action,5)==0)//Choix N°5
			{
				mpz_set_ui(exp,4423);
			}
		else if(mpz_cmp_ui(action,6)==0)//Choix N°6
			{
				mpz_set_ui(exp,9941);
			}
		else if(mpz_cmp_ui(action,7)==0)//Choix N°7
			{
				mpz_set_ui(exp,11213);
			}
		else if(mpz_cmp_ui(action,8)==0)//Choix N°8
			{
				mpz_set_ui(exp,19937);
			}
		else if(mpz_cmp_ui(action,9)==0)//Choix N°9
			{
				mpz_set_ui(exp,21701);
			}
		else if(mpz_cmp_ui(action,10)==0)//Choix N°10
			{
				mpz_set_ui(exp,23209);
			}
			gmp_printf("Vous avez choisi:(%Zd^%Zd)-1\n",nombre,exp);
			deux_pow_exp(nombre,exp);//a=2^exp
			mpz_sub_ui(nombre,nombre,1);//a=a-1 donc a=(2^exp)-1
			}
			
		else if(mpz_cmp_ui(action,11)==0)//Choix N°11
		{
			gmp_printf("**Entrer un Nombre:");
			gmp_scanf("%Zd",choix);	
			gmp_printf("Vous avez choisi:%Zd\n",choix);
			mpz_set(nombre,choix);
		}
			final=solovay_strassen(nombre,iteration);//Solovay-Strassen
			if(final==1)
			{
				gmp_printf("Test De Solovay-Strassen:Ce nombre est Probablement Premier\n");	
			}
	 
			else if(final==0)
			{
				gmp_printf("Test De Solovay-Strassen:Ce nombre est Composé\n");
			}
	
	}
	
	mpz_clears(nombre,res_jacobi,choix,ja,jn,iteration,i,i2,exp,exp2,action,NULL);//Liberation Memoire nombre,exp,iteration

	return 0;
}
