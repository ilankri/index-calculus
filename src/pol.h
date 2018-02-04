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
pol *init_pol(int deg);

// libere tout l'espace occupe par [P]
void clear_pol(pol *P);

// libere tout l'espace occupe par [l]
void clear_pol_list(pol_list *l);

// renvoie une copie du polynome [P]
pol *copie_pol(pol *P);

// renvoie 1 si [P] = 0 et 0 sinon 
int is_pol_nul(pol *P);

// renvoie 1 si [P] = [Q], 0 sinon
int egal_pol(pol *P, pol *Q);

// reduit les coefficients de [P] modulo p
void mod_pol(pol *P);

// renvoie -[P]
pol *opp_pol(pol *P);

// renvoie le polynome [coeff]X^[deg]
//pol *monome(int deg, mpz_t coeff);

// addition de 2 polynomes
pol *add_pol(pol *P, pol *Q);

// multiplication de 2 polynomes
pol *mult_pol(pol *P, pol *Q);

// renvoie [P]-[Q]
pol *soustr_pol(pol *P, pol*Q);

/*
  renvoie le reste et le quotient de la division euclidienne de [P] par [Q]; 
  si on note res la valeur de retour de cette fonction, alors res[0] designe
  le reste et res[1] le quotient  
*/ 
pol **div_eucl_pol(pol *P, pol *Q);

// affiche un polynome
void print_pol(pol *P);
