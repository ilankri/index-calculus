let main () =
  let () =
    Matrix.print @@ Matrix.make Z.[[one; of_int 2]; [of_int 3; of_int 4]];
    print_newline ()
  in
  let module P = Polynomial.Make (Q) (struct let var = "X" end) in
  (* foncteur for Fp field *)
  let module F2 = Polynomial.Make
      (struct
        include Z
        let p = of_int 2
        let modp x = rem x p
        let equal x y = equal (modp x) (modp y)
        let ( + ) x y = modp (x + y)
        let ( * ) x y = modp (x * y)
        let ( ~- ) x = modp (~- x)
        let inv x = invert x p
      end)
      (struct let var = "X" end)
  in
  let p = P.unsafe_of_list Q.[one ; one; of_int 2]
  and q = P.unsafe_of_list Q.[of_int 3; zero; of_int 4; of_int 9] in
  P.print p;
  print_newline ();
  P.print q;
  print_newline ();
  P.print @@ P.add p q;
  print_newline ();
  P.print @@ P.mul p q;
  print_newline ();
  let q, r = P.div_rem q p in
  P.print q;
  print_newline ();
  P.print r;
  print_newline ()

let () = main ()
