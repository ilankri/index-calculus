module List = struct
  (** Fill the shortest list with [None]s to equalize the length of the
      two lists.  *)
  let equalize_lengths (l : 'a list) (l' : 'b list) :
    'a option list * 'b option list =
    let rec loop (acc, acc') l l' =
      match l, l' with
      | [], [] -> List.rev acc, List.rev acc'
      | h :: t, h' :: t' -> loop (Some h :: acc, Some h' :: acc') t t'
      | h :: t, [] -> loop (Some h :: acc, None :: acc') t []
      | [], h' :: t' -> loop (None :: acc, Some h' :: acc') [] t'
    in
    loop ([], []) l l'
end
