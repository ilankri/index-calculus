#include <stdlib.h>
#include <gmp.h>
#include "pol.h"
#include "mat.h"

// polynome irreductible tel que : Fp[X]/(F) = Fq  
pol *F;

// generateur du groupe mutiplicatif Fq* 
pol *G;

// designe le degre maximal d'un polynome de la base de facteurs
int borne;

// base de polynomes irreductibles de degre inferieur ou egal a borne
pol_list *base;

// cardinal de la base de facteurs
int card_base;

// derniere puissance essayee 
mpz_t pow_prec;

// derniere puissance du generateur calculee
pol *G_pow_prec;

/*
  renvoie le polynome unitaire de degre [deg] dont les coefficients sont
  determines en decomposant [x] en base p
  Remarques : 
  1) [deg] >= 1
  2) 0 <= [x] <= p^([deg]-1)-1
*/
pol *pol_unit(int deg, mpz_t x){
  pol *P = init_pol(deg);
  mpz_t n;
  int i = 0;
  
  mpz_init_set(n, x);
  mpz_set_ui(P->coeff[deg], 1);
  while(mpz_cmp_ui(n, 0))
    mpz_tdiv_qr(n, P->coeff[i++], n, p);
  mpz_clear(n);
  return P;
}

// renvoie la liste des polynomes unitaires de degre [deg] >= 1
pol_list *pol_deg_fixe(int deg){
  pol_list *res = malloc(sizeof(pol_list));
  pol_list *pol_suiv = NULL;
  mpz_t i;
  mpz_t p_pow_deg;
  
  mpz_inits(i, p_pow_deg, NULL);
  mpz_pow_ui(p_pow_deg, p, deg);
  mpz_set(i, p_pow_deg);
  mpz_sub_ui(i, i, 1);
  while(mpz_cmp_ui(i, 0) > 0){
    pol_list *pol_courant = malloc(sizeof(pol_list));
    
    pol_courant->P = pol_unit(deg, i);
    pol_courant->pol_suiv = pol_suiv;
    pol_suiv = pol_courant;
    mpz_sub_ui(i, i, 1);
  }
  res->P = init_pol(deg);
  mpz_set_ui(res->P->coeff[deg], 1);
  res->pol_suiv = pol_suiv;
  mpz_clears(i, p_pow_deg, NULL);
  return res;
}

// renvoie tous les polynomes unitaires de degre inferieur ou egal a borne
pol_list **polynomes(){
  pol_list **res = malloc(borne*sizeof(pol_list *));
  int i;
  
  for(i=0; i<borne; i++)
    res[i] = pol_deg_fixe(i+1);
  return res;
}

/*
  renvoie un pointeur sur la liste dont la tete est le polynome
  irreductible immediatement apres [irred_courant] dans [pols]
  NB : renvoie NULL si [irred_courant] est le dernier polynome
  irreductible de [pols]
*/
pol_list *irred_suiv(pol_list **pols, pol_list *irred_courant){
  pol_list *next_pol = irred_courant->pol_suiv;
  int deg_courant = irred_courant->P->deg;
  
  if(next_pol==NULL && deg_courant==borne)
    return NULL;
  if(next_pol == NULL)
    next_pol = pols[deg_courant];
  while(next_pol->P == NULL){
    next_pol = next_pol->pol_suiv;
    if(next_pol == NULL){
      if(deg_courant == borne)
	return NULL;
      else
	next_pol = pols[deg_courant];
    }
  }
  return next_pol;
}

/* 
   met a NULL tous les polynomes de [pols] divisibles par 
   le polynome irreductible [pol_irred]
*/
void cribler(pol_list **pols, pol *pol_irred){
  int ligne_a_cribler = pol_irred->deg;
  pol_list *pol_a_cribler;
  
  while(ligne_a_cribler < borne){
    pol_a_cribler = pols[ligne_a_cribler++];
    while(pol_a_cribler != NULL){
      if(pol_a_cribler->P != NULL){
	pol **d = div_eucl_pol(pol_a_cribler->P, pol_irred);
	
	if(is_pol_nul(d[0]))
	  pol_a_cribler->P = NULL;
	clear_pol(d[0]);
	clear_pol(d[1]);
	free(d);
      }
      pol_a_cribler = pol_a_cribler->pol_suiv;
    }
  }
}

// cette fonction genere la base de facteurs
void generer_base(){
  pol_list **pols = polynomes();
  pol_list *irred_courant = pols[0];
  pol_list *next_irred;
  pol_list *facteur_suiv = NULL;
  int i;
  
  base = malloc(sizeof(pol_list));
  cribler(pols, irred_courant->P);
  while((next_irred = irred_suiv(pols, irred_courant)) != NULL){
    if(!egal_pol(irred_courant->P, G)){
      pol_list *facteur_courant = malloc(sizeof(pol_list));
      
      card_base++;
      facteur_courant->P = copie_pol(irred_courant->P);
      facteur_courant->pol_suiv = facteur_suiv;
      facteur_suiv = facteur_courant;
    }
    cribler(pols, next_irred->P);
    irred_courant = next_irred;
  }
  card_base++;
  if(!egal_pol(irred_courant->P, G)){
    base->P = copie_pol(irred_courant->P);
    base->pol_suiv = facteur_suiv;
  }
  else
    base = facteur_suiv;
  for(i=0; i<borne; i++)
    clear_pol_list(pols[i]);
  free(pols);
}

// afficher les elements composants la base de facteurs
void print_base(){
  pol_list *facteur_courant = base;

  gmp_printf("cardinal de la base : %d\n", card_base);
  while(facteur_courant != NULL){
    print_pol(facteur_courant->P);
    facteur_courant = facteur_courant->pol_suiv;
  }
}

/* 
   fonction qui factorise (dans la base de facteurs) G^[k] si 
   le second parametre est NULL, [P]G^[k] sinon
   revoie les coordonnees dans la base dans le format adapte a la partie
   algebre lineaire
   NB : cette fonction renvoie NULL s'il n'y a pas de factorisation
*/
ligne *factoriser(mpz_t k, pol *P){
  ligne *res = malloc(sizeof(ligne));
  ligne *coord_suiv = NULL;
  ligne *res_bis;
  mpz_t i;
  pol_list *facteur_courant = base;
  pol *Q; // G^k si P == NULL, PG^k sinon
  pol *R, *S;
  pol **d;
  
  if(P == NULL)
    Q = copie_pol(G_pow_prec);
  else{
    S = mult_pol(P, G_pow_prec);
    d = div_eucl_pol(S, F);
    Q = d[0];
    clear_pol(d[1]);
    free(d);
    clear_pol(S);
  }
  mpz_init_set(i, pow_prec);  
  while(mpz_cmp(i, k) < 0){
    R = copie_pol(Q);
    clear_pol(Q);
    S = mult_pol(R, G);
    d = div_eucl_pol(S, F);
    Q = d[0];
    clear_pol(d[1]);
    free(d);
    clear_pol(S);
    clear_pol(R);
    mpz_add_ui(i, i, 1);
    
  }
  mpz_set(pow_prec, k);
  clear_pol(G_pow_prec);
  G_pow_prec = copie_pol(Q);
  while(facteur_courant != NULL){
    ligne *coord_courante = malloc(sizeof(ligne));
    
    mpz_init(coord_courante->nb);
    mpz_set_ui(i, 0);
    d = div_eucl_pol(Q, facteur_courant->P);
    S = d[0];
    clear_pol(d[1]);
    free(d);
    while(is_pol_nul(S)){
      R = copie_pol(Q);
      clear_pol(Q);
      d = div_eucl_pol(R, facteur_courant->P);
      Q = d[1];
      clear_pol(d[0]);
      free(d);
      clear_pol(R);
      clear_pol(S);
      d = div_eucl_pol(Q, facteur_courant->P);
      S = d[0];
      clear_pol(d[1]);
      free(d);
      mpz_add_ui(i, i, 1);
    }
    clear_pol(S);
    mpz_set(coord_courante->nb, i);
    coord_courante->suiv = coord_suiv;
    coord_suiv = coord_courante;
    facteur_courant = facteur_courant->pol_suiv;		       
  }
  mpz_set_ui(i, 0);
  d = div_eucl_pol(Q, G);
  S = d[0];
  clear_pol(d[1]);
  free(d);
  while(is_pol_nul(S)){
    R = copie_pol(Q);
    clear_pol(Q);
    d = div_eucl_pol(R, G);
    Q = d[1];
    clear_pol(d[0]);
    free(d);
    clear_pol(R);
    clear_pol(S);
    d = div_eucl_pol(Q, G);
    S = d[0];
    clear_pol(d[1]);
    free(d);
    mpz_add_ui(i, i, 1);
  }
  clear_pol(S);
  if(Q->deg != 0){
    mpz_clear(i);
    clear_pol(Q);
    clear_ligne(coord_suiv);
    free(res);
    return NULL;
  }
  if(P == NULL){
    mpz_init_set(res->nb, k);
    mpz_sub(res->nb, res->nb, i);
  }
  else{
    mpz_init_set(res->nb, i);
    mpz_sub(res->nb, res->nb, k);
  }
  mpz_mod(res->nb, res->nb, q_moins_un); 
  res->suiv = coord_suiv;
  mpz_clear(i);
  clear_pol(Q);
  res_bis = ligne_inv(res);
  clear_ligne(res);
  return res_bis;
}
