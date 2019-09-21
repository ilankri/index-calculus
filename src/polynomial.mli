(* p designe la caracteristique de Fq  *)
module Make
    (Ring : sig
       type t
       val zero : t
       val one : t
       val equal : t -> t -> bool
       val ( + ) : t -> t -> t
       val ( * ) : t -> t -> t
       val ( ~- ) : t -> t
       val inv : t -> t
       val to_string : t -> string
     end)
    (Conf : sig val var : string end) : sig
  (* type representant un polynome.  Remarque : le polynome nul est de
     degre 0.  *)
  type t

  val zero : t

  val one : t

  val equal : t -> t -> bool

  (* renvoie -[p] *)
  val neg : t -> t

  val deg : t -> int

  val add : t -> t -> t

  val mul : t -> t -> t

  val sub : t -> t -> t

  (* renvoie le reste et le quotient de la division euclidienne de [p] par
     [q]; si on note [res] la valeur de retour de cette fonction, alors
     [fst res] designe le quotient et [snd res] le reste *)
  val div_rem : t -> t -> t * t

  val print : t -> unit

  val unsafe_of_list : Ring.t list -> t
end
