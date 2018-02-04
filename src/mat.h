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

void clear_ligne(ligne *l);

// renvoie la ligne inverse de [l], par exemple 1,2,3 devient 3,2,1
ligne *ligne_inv(ligne *l);

// renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon
int is_ligne_nulle(ligne *l);

// renvoie 1 si les lignes [l1] et [l2] sont egales, 0 sinon
int egal_ligne(ligne *l1, ligne *l2);

// affiche une ligne
void print_ligne(ligne *l);

// affiche une matrice
void print_mat(mat *M);
