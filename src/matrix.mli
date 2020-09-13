type t

val make : Z.t list list -> t

val print : t -> unit

val nth_line : int -> t -> Vector.t
(* les lignes sont numerotés de 0 a nb_ligne-1 *)
(* peme_ligne(p,M) renvoie la p+1eme ligne de M *)
