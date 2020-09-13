type t = Z.t list

let make l = l

let rev = List.rev

(* renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon *)
let is_zero v = List.for_all (Z.(equal zero)) v

let equal v v' =
  try List.for_all2 (fun x x' -> Z.equal x x') v v' with
  | Invalid_argument _ -> false

let print v = print_string @@ String.concat " " (List.map Z.to_string v)

let scal_mult k v = List.map (Z.mul k) v

(* TODO: Handle the case where vector are not of the same length.  *)
let add v v' = List.map2 Z.add v v'

let set_nth n i v = ExtStdlib.List.replace_nth n i v
