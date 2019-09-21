type t = Z.t list

let make l = l

let rev = List.rev

(* renvoie 1 si la ligne [l] est une ligne de zeros, 0 sinon *)
let is_zero v = List.for_all (Z.(equal zero)) v

let equal v v' =
  try List.for_all2 (fun x x' -> Z.equal x x') v v' with
  | Invalid_argument _ -> false

let print v = print_string @@ String.concat " " (List.map Z.to_string v)
