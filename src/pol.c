#include <stdlib.h>
#include <gmp.h>

// p designe la caracteristique de Fq
mpz_t p;

/*
   type representant un polynome
   Remarques :
   1) pour tout 0 <= i <= deg, coeff[i] designe le coefficient devant X^i
   2) le polynome nul est de degre 0
*/
typedef struct{
  int deg;
  mpz_t *coeff;
} pol;

// type liste de polynomes
typedef struct pol_list{
  pol *P;
  struct pol_list *pol_suiv;
} pol_list;

// renvoie un polynome de degre [deg] en mettant tous ses coefficients a 0
pol *init_pol(int deg){
  int i;
  pol *P = malloc(sizeof(pol));

  P->deg = deg;
  P->coeff = malloc((deg+1)*sizeof(mpz_t));
  for(i=0; i<=deg; i++)
    mpz_init(P->coeff[i]);
  return P;
}

// libere tout l'espace occupe par [P]
void clear_pol(pol *P){
  if(P != NULL){
    int i;

    for(i=0; i<=P->deg; i++)
      mpz_clear(P->coeff[i]);
    free(P->coeff);
    free(P);
  }
}

// libere tout l'espace occupe par [l]
void clear_pol_list(pol_list *l){
  if(l != NULL){
    clear_pol(l->P);
    clear_pol_list(l->pol_suiv);
    free(l);
  }
}

// renvoie une copie du polynome [P]
pol *copie_pol(pol *P){
  pol *Q = init_pol(P->deg);
  int i;

  for(i=0; i<=Q->deg; i++)
    mpz_set(Q->coeff[i], P->coeff[i]);
  return Q;
}

// renvoie 1 si [P] = 0 et 0 sinon
int is_pol_nul(pol *P){
  if(P->deg > 0)
    return 0;
  return mpz_cmp_ui(P->coeff[0], 0) == 0;
}

// renvoie 1 si [P] = [Q], 0 sinon
int egal_pol(pol *P, pol *Q){
  int i;

  if(P->deg != Q->deg)
    return 0;
  for(i=0; i<=P->deg; i++)
    if(mpz_cmp(P->coeff[i], Q->coeff[i]))
      return 0;
  return 1;
}

// reduit les coefficients de [P] modulo p
void mod_pol(pol *P){
  int i;

  for(i=0; i<=P->deg; i++)
    mpz_mod(P->coeff[i], P->coeff[i], p);
}

// renvoie -[P]
pol *opp_pol(pol *P){
  pol *Q = init_pol(P->deg);
  int i;

  for(i=0; i<=Q->deg; i++)
    mpz_neg(Q->coeff[i], P->coeff[i]);
  mod_pol(Q);
  return Q;
}

// addition de 2 polynomes
pol *add_pol(pol *P, pol *Q){
  pol *R;
  int i;

  if(P->deg == Q->deg){
    pol *S;

    R = init_pol(P->deg);
    for(i=0; i<=P->deg; i++)
      mpz_add(R->coeff[i], P->coeff[i], Q->coeff[i]);
    mod_pol(R);
    while(R->deg>=1 && (mpz_cmp_ui(R->coeff[R->deg], 0) == 0))
      R->deg--;
    S = copie_pol(R);
    clear_pol(R);
    return S;
  }
  else{
    int deg_min;

    if(P->deg > Q->deg){
      deg_min = Q->deg;
      R = init_pol(P->deg);
      for(i=deg_min+1; i<=P->deg; i++)
	mpz_set(R->coeff[i], P->coeff[i]);
    }
    else{
      deg_min = P->deg;
      R = init_pol(Q->deg);
      for(i=deg_min+1; i<=Q->deg; i++)
	mpz_set(R->coeff[i], Q->coeff[i]);
    }
    for(i=0; i<=deg_min; i++){
      mpz_add(R->coeff[i], P->coeff[i], Q->coeff[i]);
      mpz_mod(R->coeff[i], R->coeff[i], p);
    }
    return R;
  }
}

// multiplication de 2 polynomes
pol *mult_pol(pol *P, pol *Q){
  if(is_pol_nul(P) ||  is_pol_nul(Q)){
    return init_pol(0);
  }
  else{
    pol *R = init_pol(P->deg+Q->deg);;
    int i, j;

    for(i=0; i<=P->deg; i++)
      for(j=0; j<=Q->deg; j++)
	mpz_addmul(R->coeff[i+j], P->coeff[i], Q->coeff[j]);
    mod_pol(R);
    return R;
  }
}

// renvoie [P]-[Q]
pol *soustr_pol(pol *P, pol*Q){
  return add_pol(P, opp_pol(Q));
}

/*
  renvoie le reste et le quotient de la division euclidienne de [P] par [Q];
  si on note res la valeur de retour de cette fonction, alors res[0] designe
  le reste et res[1] le quotient
*/
pol **div_eucl_pol(pol *P, pol *Q){
  pol **res = malloc(2*sizeof(pol *));

  res[0] = copie_pol(P);
  if(P->deg < Q->deg){
    res[1] = init_pol(0);
    return res;
  }
  res[1] = init_pol(P->deg-Q->deg);
  while(res[0]->deg >= Q->deg){
    pol *S, *T;
    int deg_monome = res[0]->deg - Q->deg;
    mpz_t coeff_monome;
    mpz_t inv_coeff_Q; // l'inverse du coefficient dominant de [Q]
    int i;

    mpz_init(coeff_monome);
    mpz_init(inv_coeff_Q);
    mpz_invert(inv_coeff_Q, Q->coeff[Q->deg], p);
    mpz_mul(coeff_monome, P->coeff[P->deg], inv_coeff_Q);
    mpz_mod(coeff_monome, coeff_monome, p);
    mpz_set(res[1]->coeff[deg_monome], coeff_monome);
    S = init_pol(Q->deg + deg_monome);
    for(i=0; i<=Q->deg; i++){
      mpz_mul(S->coeff[i+deg_monome], Q->coeff[i], coeff_monome);
      mpz_mod(S->coeff[i+deg_monome], S->coeff[i+deg_monome], p);
    }
    T = copie_pol(res[0]);
    clear_pol(res[0]);
    res[0] = soustr_pol(T, S);
    mpz_clears(coeff_monome, inv_coeff_Q, NULL);
    clear_pol(S);
    clear_pol(T);
  }
  return res;
}

// affiche un polynome
void print_pol(pol *P){
  int i;

  for(i=P->deg; i>=0; i--)
    gmp_printf("a%d=%Zd ", i, P->coeff[i]);
  gmp_printf("\n");
}
