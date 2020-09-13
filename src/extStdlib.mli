module List : sig
  val equalize_lengths : 'a list -> 'b list -> 'a option list * 'b option list
  (** Fill the shortest list with [None]s to equalize the length of the
      two lists.  *)

  val replace_nth : int -> 'a -> 'a list -> 'a list
end
