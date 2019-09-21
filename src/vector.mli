type t

val make : Z.t list -> t

val rev : t -> t

(* renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon *)
val is_zero : t ->  bool

val equal : t -> t -> bool

val print : t -> unit
