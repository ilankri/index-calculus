module List = struct
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

  let rec replace_nth n x = function
    | [] -> []
    | h :: t -> if n = 0 then x :: t else h :: replace_nth (n - 1) x t
end
