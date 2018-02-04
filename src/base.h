// polynome irreductible tel que : Fp[X]/(F) = Fq  
pol *F;

// generateur du groupe mutiplicatif Fq* 
pol *G;

// l'ordre de Fq*
mpz_t q_moins_un;

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

// cette fonction genere la base de facteurs
void generer_base();

// afficher les elements composants la base de facteurs
void print_base();

/* 
   fonction qui factorise (dans la base de facteurs) G^[k] si 
   le second parametre est NULL, [P]G^[k] sinon
   revoie les coordonnees dans la base dans le format adapte a la partie
   algebre lineaire
   NB : cette fonction renvoie NULL s'il n'y a pas de factorisation
*/
ligne *factoriser(mpz_t k, pol *P);
