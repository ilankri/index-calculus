// programme principal
// compilation : gcc pol.c mat.c base.c gauss.c calc_index.c -lgmp

#include <stdlib.h>
#include <gmp.h>
#include "pol.h"
#include "mat.h"
#include "base.h"
#include "gauss.h"

// dimension de Fq vu comme Fp-espace vectoriel
int n;

/* 
   designe le nombre de relations a generer entre les logarithmes 
   discrets des elements de la base
*/
int nb_rel;

// liste des logarithmes discrets des elements de la base
ligne *log_base;

// element de Fq* dont on cherche le logarithme discret 
pol *H; 

// matrice representant le systeme lineaire a resoudre
mat *S;

// le logarithme de H 
mpz_t log_H; 

/*
  cette fonction genere le systeme d'equations lineaires dont les 
  inconnues sont les logarithmes discrets des elements de la base
*/
mat *generer_syst(){
  mat *M = NULL;
  int i = nb_rel;
  ligne *l = NULL;
  ligne *l_prec = NULL;
  mpz_t premier_suiv;
  
  mpz_init(premier_suiv);
  while(i--){
    mat *N = malloc(sizeof(mat));
    
    N->nb_ligne = i;
    l = NULL;
    while(l == NULL || is_ligne_nulle(l) || egal_ligne(l, l_prec)){
      clear_ligne(l);
      mpz_nextprime(premier_suiv, premier_suiv);
      l = factoriser(premier_suiv, NULL);
    }
    N->l = l;
    l_prec = l;
    N->suiv = M;
    M = N;
  }
  mpz_clear(premier_suiv);
  return M;
}

/* 
   cette fonction trouve s tel que (G^s)H se decompose dans la base de 
   facteur et en deduit le log discret de H
*/
void log_discret(){
  ligne *l = NULL;
  mpz_t premier_suiv;
  ligne *log_base_suiv = log_base; 
  
  mpz_init(premier_suiv);
  mpz_init(log_H);
  l = factoriser(premier_suiv, H);
  if(l == NULL || is_ligne_nulle(l)){
    mpz_set_ui(premier_suiv, 1);
    l = factoriser(premier_suiv, H);
  }
  while(l == NULL || is_ligne_nulle(l)){
    clear_ligne(l);
    mpz_nextprime(premier_suiv, premier_suiv);
    l = factoriser(premier_suiv, H);
  }
  while(log_base_suiv != NULL){
    mpz_addmul(log_H, log_base_suiv->nb, l->nb);
    mpz_mod(log_H, log_H, q_moins_un);
    log_base_suiv = log_base_suiv->suiv;
    l = l->suiv;
  }
  mpz_add(log_H, log_H, l->nb);
  mpz_mod(log_H, log_H, q_moins_un);
  
}

/* Test dans F_2^17 = < X >, on cherche le log discret de 
   H = X^10 + X^5 + X^2 */
int main(int argc, char *argv[]){
  // initialisation des parametres
  n = 17;
  mpz_inits(p, q_moins_un, pow_prec, NULL);
  mpz_set_ui(p, 2);
  mpz_pow_ui(q_moins_un, p, n);
  mpz_sub_ui(q_moins_un, q_moins_un, 1);
  F = init_pol(n);
  mpz_set_ui(F->coeff[0], 1);
  mpz_set_ui(F->coeff[14], 1);
  mpz_set_ui(F->coeff[n], 1);
  G = init_pol(1);
  mpz_set_ui(G->coeff[1], 1);
  H = init_pol(10);
  mpz_set_ui(H->coeff[2], 1);
  mpz_set_ui(H->coeff[5], 1);
  mpz_set_ui(H->coeff[10], 1);
  mpz_set_ui(pow_prec, 1);
  G_pow_prec = copie_pol(G);
  borne = 3;
  
  // construction de la base
  generer_base();
  
  // generation du systeme lineaire
  nb_rel = 2*card_base;
  S = generer_syst();
  
  // resolution du systeme
  m = card_base + 1;
  log_base = piv_gauss(S, nb_rel);
  
  // calcul final
  mpz_set_ui(pow_prec, 0);
  G_pow_prec = init_pol(0);
  mpz_set_ui(G_pow_prec->coeff[0], 1);
  log_discret();
  gmp_printf("F(X) : ");
  print_pol(F);
  gmp_printf("G(X) : ");
  print_pol(G);
  gmp_printf("F_%Zd^%d = F_%Zd[X]/(F(X)) = < G(X) >\n\n", p, n, p);
  gmp_printf("H(X) : ");
  print_pol(H);
  gmp_printf("log(H(X)) = %Zd\n", log_H);
  return 0;
}
