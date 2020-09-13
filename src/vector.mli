type t

val make : Z.t list -> t

val rev : t -> t

(* renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon *)
val is_zero : t ->  bool

val equal : t -> t -> bool

val print : t -> unit

val scal_mult : Z.t -> t -> t

val add : t -> t -> t

val set_nth : int -> Z.t -> t -> t
(* les éléments d'une ligne sont numérotés de 0 a m-1 *)
(* peme_elem(p,n,l) met le p+1eme element de l dans n *)
