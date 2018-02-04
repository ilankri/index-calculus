#include <stdlib.h>
#include <gmp.h>

// l'ordre de Fq*
mpz_t q_moins_un;

typedef struct ligne{
  mpz_t nb;
  struct ligne *suiv;
} ligne;

typedef struct mat{
  ligne *l;
  struct mat *suiv;
  int nb_ligne;
} mat;

// libere tout l'espace occupe par [l]
void clear_ligne(ligne *l){
  if(l != NULL){
    mpz_clear(l->nb);
    clear_ligne(l->suiv);
    free(l);
  }
}

// renvoie la ligne inverse de [l], par exemple 1,2,3 devient 3,2,1
ligne *ligne_inv(ligne *l){
  ligne *res = NULL;
  ligne *l_next = l;
  
  while(l_next != NULL){
    ligne *res_courant = malloc(sizeof(ligne));
    
    res_courant->suiv = res;
    mpz_init_set(res_courant->nb, l_next->nb);
    res = res_courant;
    l_next = l_next->suiv;
  }
  return res;
}

// renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon
int is_ligne_nulle(ligne *l){
  ligne *l_bis = l;
  
  if(l == NULL)
    return 0;
  while(l_bis != NULL){
    if(mpz_cmp_ui(l_bis->nb, 0))
      return 0;
    l_bis = l_bis->suiv;
  }
  return 1;
}

/* 
   renvoie 1 si les lignes (de meme longueur) [l1] et [l2] sont egales, 
   0 sinon
*/
int egal_ligne(ligne *l1, ligne *l2){
  ligne *l1_bis = l1;
  ligne *l2_bis = l2;
  
  if(l1==NULL || l2==NULL)
    return 0;
  while(l1_bis != NULL){
    if(mpz_cmp(l1_bis->nb, l2_bis->nb))
      return 0;
    l1_bis = l1_bis->suiv;
    l2_bis = l2_bis->suiv;
  }
  return 1;
}

// affiche une ligne
void print_ligne(ligne *l){
  ligne *l_bis = l;
  
  while(l_bis != NULL){
    gmp_printf("%Zd ", l_bis->nb);
    l_bis = l_bis->suiv;
  }
}

// affiche une matrice
void print_mat(mat *M){
  mat *N = M;
  
  while(N != NULL){
    gmp_printf("%d   ", N->nb_ligne);
    print_ligne(N->l);
    gmp_printf("\n");
    N = N->suiv;
  }  
}
