#include <stdlib.h>
#include <gmp.h>
#include "mat.h"
#include "pol.h"
#include "base.h"

int m;    //nb de collones +1 :le second membre

void perm_circ(int k, mat *M, int nb_lignes_M) {
  int i;
  ligne* L=peme_ligne(k, M);
  mat *A=M;
  for(i=0;i<k;i++)
    A=A->suiv;
  for(i=k;i<nb_lignes_M-1;i++) {
    A->l=(A->suiv)->l;
    A=A->suiv;
  }
  A->l=L;
}

ligne *piv_gauss(mat *M, int nb_lignes_M) {
  int i,j,k;
  ligne *ieme;
  ligne *jeme;
  mpz_t N, P;
  int c = 0;
  mpz_inits(N, P, NULL);
  for(i=0;i<m-2;i++) {
    ieme=peme_ligne(i, M);
    peme_elem(i, N, ieme);
    while(mpz_cmp_ui(N, 0)==0 && c++<nb_lignes_M-1-i) {
      perm_circ(i,M,nb_lignes_M);
      ieme=peme_ligne(i,M);
      peme_elem(i, N, ieme);
    }
    peme_elem(i, P, ieme);
    for(j=i+1;j<nb_lignes_M;j++) {
      jeme=peme_ligne(j, M);
      peme_elem(i, N,  jeme);
      if(mpz_cmp_ui(N,0)!=0) {
	mpz_invert(N, N, q_moins_un);
	mpz_mul(N, N, P);
	mpz_neg(N, N);
	mult_scal(N, jeme);
	add_lignes(jeme, ieme);
      }
    }
    jeme=peme_ligne(i+1,M);
    peme_elem(i+1, N, jeme);
    c = 0;
    //c pour arret (pas infini)
    while(mpz_cmp_ui(N, 0)==0 && c++<nb_lignes_M-2-i) {
      perm_circ(i+1, M, nb_lignes_M);
      jeme=peme_ligne(i+1,M);
      peme_elem(i+1, N, jeme);
    }
  }
  ligne *L;
  ligne *R;
  mpz_t I,J,K;
  mpz_init(I);
  mpz_init(J);
  mpz_init(K);
  L = alloue_ligne();
  for(i=m-2;i>=0;i--) {
    ieme=peme_ligne(i, M);
    peme_elem(m-1, N, ieme);
    mpz_neg(N, N);
    for(j=i+1;j<m-1;j++) {
      peme_elem(j, I, L);
      peme_elem(j, J, ieme);
      mpz_addmul(N, I, J);
    }
    mpz_neg(N, N);
    peme_elem(i, I, ieme);
    mpz_set(K, I);
    mpz_invert(K, K, q_moins_un);
    peme_elem(i, I, L);
    R=L;
    for(k=0;k<i;k++)
      R=R->suiv;
    mpz_mul(N, K, N);
    mpz_mod(R->nb, N, q_moins_un);
  }
  return L;
}
