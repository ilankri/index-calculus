let main () =
  let () =
    Matrix.print @@ Matrix.make Z.[[one; of_int 2]; [of_int 3; of_int 4]];
    print_newline ()
  in
  (* default module for default conf *)
  let module P = Polynomial.Make (Q) (struct let var = "X" end) in
  (* foncteur for Fp field *)
  let module PP = Polynomial.Make
      (struct
        include Z
        let p = of_int 2
        let modp x = rem x p
        let equal x y = equal (modp x) (modp y)
        let add x y = modp @@ add x y
        let mul x y = modp @@ mul x y
        let neg x =
          let res = modp @@ neg x in
          if sign res = -1 then res + p else res
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
  print_newline ();
  let p = PP.unsafe_of_list Z.[one ; one]
  and q = PP.unsafe_of_list Z.[one; zero; zero; one] in
  PP.print p;
  print_newline ();
  PP.print q;
  print_newline ();
  PP.print @@ PP.add p q;
  print_newline ();
  PP.print @@ PP.mul p q;
  print_newline ();
  let q, r = PP.div_rem q p in
  PP.print q;
  print_newline ();
  PP.print r

let () = main ()
