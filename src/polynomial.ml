(* It is not a ring, it is a field because of inv! Two functors one for
   polynomials on a field and another one polynomials on a ring.*)

module Make
    (Ring : sig
       type t
       val zero : t
       val one : t
       val equal : t -> t -> bool
       val add : t -> t -> t
       val mul : t -> t -> t
       val neg : t -> t
       val inv : t -> t
       val to_string : t -> string
     end)
    (Conf : sig val var : string end) = struct
  (* type representant un polynome.  Remarques :

     1) pour tout 0 <= i <= deg, coeff[i] designe le coefficient devant
        X^i

     2) le polynome nul est de degre 0 *)
  type t = Ring.t list

  let zero = [Ring.zero]

  let one = [Ring.one]

  let equal p q =
    try List.for_all2 Ring.equal p q with Invalid_argument _ -> false

  let neg p = List.map Ring.neg p

  (* deg(0) = 0 or not defined? deg(k) = 1 where k <> 0?  *)
  let deg p = List.length p - 1

  let add p q =
    let r =
      let p, q = ExtStdlib.List.equalize_lengths p q in
      try
        List.map2 (fun a b ->
          Ring.add
            (BatOption.default Ring.zero a) (BatOption.default Ring.zero b)
        ) p q
      with
      | Invalid_argument _ -> assert false
    in
    match
      List.rev @@ BatList.drop_while (Ring.equal Ring.zero) (List.rev r)
    with
    | [] -> zero
    | _ :: _ as r -> r

  let sub p q = add p (neg q)

  let monomial (i : int) (a : Ring.t) : t = BatList.make i Ring.zero @ [a]

  let unsafe_mul_by_monomial (i : int) (a : Ring.t) (p : t) : t =
    BatList.make i Ring.zero @ List.map (Ring.mul a) p

  let mul p q =
    BatList.fold_lefti (fun acc i a ->
      add acc (unsafe_mul_by_monomial i a q)
    ) zero p

  let get_coeff (i : int) (p : t) : Ring.t = List.nth p i

  let div_rem p q =
    let rec loop acc p =
      let deg_p = deg p
      and deg_q = deg q in
      if deg_p < deg_q
      then (acc, p)
      else
        let monomial =
          let a =
            Ring.mul (Ring.inv @@ get_coeff deg_q q) (get_coeff deg_p p)
          and i = deg_p - deg q in
          monomial i a
        in
        loop (add acc monomial) (sub p (mul monomial q))
    in
    loop zero p

  let unsafe_of_list l = l

  let print p =
    let monomials =
      let string_of_monomial i coeff =
        let i =
          match i with
          | 0 -> ""
          | 1 -> Conf.var
          | i -> Format.sprintf "%s^%d" Conf.var i
        and coeff =
          if coeff = Ring.zero
          then None
          else if coeff = Ring.one && i <> 0
          then Some ""
          else Some (Ring.to_string coeff)
        in
        BatOption.map (fun coeff -> Format.sprintf "%s%s" coeff i) coeff
      in
      BatList.filteri_map string_of_monomial p
    in
    print_string @@ String.concat " + " (List.rev monomials)
end
