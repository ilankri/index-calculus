type t = Vector.t list

let make l = List.map Vector.make l

let print m = List.iter (fun v -> Vector.print v; print_newline ()) m
