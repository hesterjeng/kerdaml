let () = Random.self_init ()

(* Box-Muller transform: two independent standard normals *)
let box_muller () =
  let u1 = Random.float 1.0 in
  let u2 = Random.float 1.0 in
  let r = sqrt (-2.0 *. log u1) in
  let theta = 2.0 *. Float.pi *. u2 in
  (r *. cos theta, r *. sin theta)

let normal ~mu ~sigma =
  let (z, _) = box_muller () in
  mu +. sigma *. z

let () =
  (* Generate mixture of two normals: 0.4 * N(-2, 0.8) + 0.6 * N(2, 1.2) *)
  let n = 1000 in
  let data = Kerdaml.create_vec n in
  for i = 0 to n - 1 do
    let x =
      if Random.float 1.0 < 0.4 then normal ~mu:(-2.0) ~sigma:0.8
      else normal ~mu:2.0 ~sigma:1.2
    in
    data.{i} <- x
  done;

  let h = Kerdaml.compute_bandwidth Scott data in
  Printf.printf "# bandwidth (Scott): %.6f\n" h;

  let result = Kerdaml.estimate ~bw:Scott ~n_points:512 data in

  (* Print grid *)
  let np = Bigarray.Array1.dim result.points in
  Printf.printf "# x\tdensity\n";
  for i = 0 to np - 1 do
    Printf.printf "%.6f\t%.6f\n" result.points.{i} result.density.{i}
  done;

  (* Trapezoidal integral as sanity check *)
  let integral = ref 0.0 in
  for i = 0 to np - 2 do
    let dx = result.points.{i + 1} -. result.points.{i} in
    integral := !integral +. 0.5 *. (result.density.{i} +. result.density.{i + 1}) *. dx
  done;
  Printf.printf "# integral (trapezoidal): %.6f\n" !integral
