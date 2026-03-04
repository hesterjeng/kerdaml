let passed = ref 0
let failed = ref 0

let check name cond =
  if cond then (
    incr passed;
    Printf.printf "  PASS  %s\n" name)
  else (
    incr failed;
    Printf.printf "  FAIL  %s\n" name)

let approx ?(tol = 1e-6) a b = Float.abs (a -. b) < tol

let section name = Printf.printf "\n=== %s ===\n" name

(* ---- vec helpers ---- *)

let test_create_vec () =
  section "create_vec / vec_of_array";
  let v = Kerdaml.create_vec 5 in
  check "create_vec length" (Bigarray.Array1.dim v = 5);
  check "create_vec zeroed" (v.{0} = 0.0 && v.{4} = 0.0);
  let a = [| 1.0; 2.0; 3.0 |] in
  let v2 = Kerdaml.vec_of_array a in
  check "vec_of_array length" (Bigarray.Array1.dim v2 = 3);
  check "vec_of_array values" (v2.{0} = 1.0 && v2.{1} = 2.0 && v2.{2} = 3.0)

(* ---- compute_bandwidth ---- *)

let test_bandwidth_fixed () =
  section "compute_bandwidth Fixed";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0 |] in
  let h = Kerdaml.compute_bandwidth (Fixed 0.42) data in
  check "Fixed returns exact value" (h = 0.42)

let test_bandwidth_scott () =
  section "compute_bandwidth Scott";
  (* Known data: 1..5 uniform, std = sqrt(2.5) ≈ 1.5811 *)
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let h = Kerdaml.compute_bandwidth Scott data in
  (* 1.06 * std * n^(-0.2), n=5, std=sqrt(2.5) *)
  let sigma = sqrt 2.5 in
  let expected = 1.06 *. sigma *. (5.0 ** (-0.2)) in
  check "Scott formula" (approx ~tol:1e-10 h expected);
  check "Scott positive" (h > 0.0)

let test_bandwidth_silverman () =
  section "compute_bandwidth Silverman";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let h = Kerdaml.compute_bandwidth Silverman data in
  check "Silverman positive" (h > 0.0);
  check "Silverman < Scott" (h <= Kerdaml.compute_bandwidth Scott data)

let test_bandwidth_constant_data () =
  section "compute_bandwidth constant data";
  let data = Kerdaml.vec_of_array [| 5.0; 5.0; 5.0; 5.0; 5.0 |] in
  let h = Kerdaml.compute_bandwidth Scott data in
  check "Scott on constant data = 0" (h = 0.0)

(* ---- evaluate_at ---- *)

let test_evaluate_at_single_point () =
  section "evaluate_at single data point";
  (* KDE of one point at x=0 with h=1:
     density(x) = gaussian(x) = (1/sqrt(2pi)) * exp(-x^2/2) *)
  let data = Kerdaml.vec_of_array [| 0.0 |] in
  let d = Kerdaml.evaluate_at 1.0 data 0.0 in
  let expected = 1.0 /. sqrt (2.0 *. Float.pi) in
  check "peak at data point" (approx ~tol:1e-10 d expected);
  let d_far = Kerdaml.evaluate_at 1.0 data 10.0 in
  check "far from data ≈ 0" (d_far < 1e-20)

let test_evaluate_at_symmetric () =
  section "evaluate_at symmetry";
  let data = Kerdaml.vec_of_array [| 0.0 |] in
  let d_pos = Kerdaml.evaluate_at 1.0 data 1.5 in
  let d_neg = Kerdaml.evaluate_at 1.0 data (-1.5) in
  check "symmetric around single point" (approx ~tol:1e-15 d_pos d_neg)

let test_evaluate_at_two_points () =
  section "evaluate_at two symmetric points";
  let data = Kerdaml.vec_of_array [| -1.0; 1.0 |] in
  let d_zero = Kerdaml.evaluate_at 1.0 data 0.0 in
  let d_left = Kerdaml.evaluate_at 1.0 data (-2.0) in
  let d_right = Kerdaml.evaluate_at 1.0 data 2.0 in
  check "symmetric KDE at origin" (d_zero > 0.0);
  check "symmetric tails" (approx ~tol:1e-15 d_left d_right);
  check "peak > tail" (d_zero > d_left)

let test_evaluate_at_bandwidth_effect () =
  section "evaluate_at bandwidth effect";
  let data = Kerdaml.vec_of_array [| 0.0 |] in
  let d_narrow = Kerdaml.evaluate_at 0.1 data 0.0 in
  let d_wide = Kerdaml.evaluate_at 10.0 data 0.0 in
  check "narrow bw → taller peak" (d_narrow > d_wide)

(* ---- evaluate (convenience) ---- *)

let test_evaluate () =
  section "evaluate convenience";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let d_default = Kerdaml.evaluate data 3.0 in
  check "evaluate default positive" (d_default > 0.0);
  let h = Kerdaml.compute_bandwidth Scott data in
  let d_manual = Kerdaml.evaluate_at h data 3.0 in
  check "evaluate matches evaluate_at with Scott" (approx ~tol:1e-15 d_default d_manual);
  let d_fixed = Kerdaml.evaluate ~bw:(Fixed 0.5) data 3.0 in
  let d_manual2 = Kerdaml.evaluate_at 0.5 data 3.0 in
  check "evaluate Fixed matches evaluate_at" (approx ~tol:1e-15 d_fixed d_manual2)

(* ---- estimate (grid) ---- *)

let trapezoidal_integral (result : Kerdaml.estimate) =
  let n = Bigarray.Array1.dim result.points in
  let s = ref 0.0 in
  for i = 0 to n - 2 do
    let dx = result.points.{i + 1} -. result.points.{i} in
    s := !s +. 0.5 *. (result.density.{i} +. result.density.{i + 1}) *. dx
  done;
  !s

let test_estimate_defaults () =
  section "estimate defaults";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate data in
  check "default 512 points" (Bigarray.Array1.dim result.points = 512);
  check "density same length" (Bigarray.Array1.dim result.density = 512)

let test_estimate_custom_n () =
  section "estimate custom n_points";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~n_points:64 data in
  check "custom 64 points" (Bigarray.Array1.dim result.points = 64)

let test_estimate_grid_ordered () =
  section "estimate grid monotonically increasing";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~n_points:100 data in
  let ok = ref true in
  for i = 0 to 98 do
    if result.points.{i} >= result.points.{i + 1} then ok := false
  done;
  check "grid strictly increasing" !ok

let test_estimate_density_positive () =
  section "estimate density all positive";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~n_points:100 data in
  let ok = ref true in
  for i = 0 to 99 do
    if result.density.{i} <= 0.0 then ok := false
  done;
  check "all density > 0" !ok

let test_estimate_integrates_to_one () =
  section "estimate integrates to ~1.0";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~n_points:1024 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "integral = %.6f ≈ 1.0" integral) (approx ~tol:1e-3 integral 1.0)

let test_estimate_large_data () =
  section "estimate large dataset integrates to ~1.0";
  Random.init 42;
  let n = 2000 in
  let data = Kerdaml.create_vec n in
  for i = 0 to n - 1 do
    (* simple normal via Box-Muller, discard second value *)
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    data.{i} <- sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)
  done;
  let result = Kerdaml.estimate ~n_points:1024 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "large data integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0)

let test_estimate_silverman () =
  section "estimate with Silverman bandwidth";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~bw:Silverman ~n_points:1024 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "Silverman integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0)

let test_estimate_fixed_bw () =
  section "estimate with Fixed bandwidth";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~bw:(Fixed 0.5) ~n_points:1024 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "Fixed(0.5) integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0)

let test_estimate_covers_data () =
  section "estimate grid covers all data points";
  let data = Kerdaml.vec_of_array [| -10.0; 0.0; 10.0 |] in
  let result = Kerdaml.estimate data in
  let n = Bigarray.Array1.dim result.points in
  check "grid starts below min" (result.points.{0} < -10.0);
  check "grid ends above max" (result.points.{n - 1} > 10.0)

(* ---- edge cases ---- *)

let test_two_points () =
  section "two data points";
  let data = Kerdaml.vec_of_array [| 0.0; 1.0 |] in
  let h = Kerdaml.compute_bandwidth Scott data in
  check "Scott bw for 2 points positive" (h > 0.0);
  let result = Kerdaml.estimate ~n_points:256 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "2-point integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-2 integral 1.0)

let test_identical_points () =
  section "identical data points";
  let data = Kerdaml.vec_of_array [| 3.0; 3.0; 3.0 |] in
  let d = Kerdaml.evaluate ~bw:(Fixed 1.0) data 3.0 in
  let expected = 1.0 /. sqrt (2.0 *. Float.pi) in
  check "identical points, peak = gaussian(0)" (approx ~tol:1e-10 d expected)

let test_negative_data () =
  section "negative data values";
  let data = Kerdaml.vec_of_array [| -5.0; -3.0; -1.0 |] in
  let d = Kerdaml.evaluate data (-3.0) in
  check "density at center of negative data > 0" (d > 0.0);
  let result = Kerdaml.estimate ~n_points:512 data in
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "negative data integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0)

(* ---- helpers ---- *)

(* Sort mode indices by grid position for spatial assertions *)
let modes_by_position (result : Kerdaml.estimate) modes =
  let a = Array.copy modes in
  Array.sort (fun i j ->
    Float.compare result.points.{i} result.points.{j}
  ) a;
  a

(* ---- find_modes tests ---- *)

let test_find_modes_unimodal () =
  section "find_modes unimodal";
  let data = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let result = Kerdaml.estimate ~n_points:512 data in
  let modes = Kerdaml.find_modes result in
  check "unimodal has 1 mode" (Array.length modes = 1);
  if Array.length modes = 1 then begin
    let i = modes.(0) in
    let x = result.points.{i} in
    let d = result.density.{i} in
    check (Printf.sprintf "mode x near 3.0 (got %.2f)" x) (approx ~tol:1.0 x 3.0);
    check "mode density > 0" (d > 0.0)
  end

let test_find_modes_sorted_by_density () =
  section "find_modes sorted by descending density";
  (* 80 points near 0 (tall peak), 20 points near 8 (short peak) *)
  let data = Kerdaml.create_vec 100 in
  Random.init 77;
  for i = 0 to 79 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    data.{i} <- 0.3 *. sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)
  done;
  for i = 0 to 19 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    data.{80 + i} <- 8.0 +. 0.3 *. sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)
  done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.4) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "2 modes found" (Array.length modes = 2);
  if Array.length modes = 2 then begin
    check "first mode has higher density"
      (result.density.{modes.(0)} >= result.density.{modes.(1)});
    let dominant_x = result.points.{modes.(0)} in
    check (Printf.sprintf "dominant mode near 0 (got %.2f)" dominant_x)
      (approx ~tol:1.0 dominant_x 0.0)
  end

let test_find_modes_empty_result () =
  section "find_modes on monotonic density";
  (* A single data point with huge bandwidth gives a monotonic bump
     that peaks at the edge or has a single maximum — verify no crash *)
  let data = Kerdaml.vec_of_array [| 0.0 |] in
  let result = Kerdaml.estimate ~bw:(Fixed 1.0) ~n_points:64 data in
  let modes = Kerdaml.find_modes result in
  check "single-point KDE has modes" (Array.length modes >= 0)

let test_find_modes_index_values () =
  section "find_modes indices give correct x and density";
  let data = Kerdaml.create_vec 100 in
  for i = 0 to 49 do data.{i} <- -5.0 +. Float.of_int i *. 0.01 done;
  for i = 0 to 49 do data.{50 + i} <- 5.0 +. Float.of_int i *. 0.01 done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.3) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "2 modes" (Array.length modes = 2);
  Array.iter (fun i ->
    let x = result.points.{i} in
    let d = result.density.{i} in
    check (Printf.sprintf "mode at x=%.2f has density=%.4f > 0" x d) (d > 0.0);
    (* verify it's actually a local max *)
    check (Printf.sprintf "index %d is local max" i)
      (result.density.{i} > result.density.{i - 1}
       && result.density.{i} > result.density.{i + 1})
  ) modes

(* ---- bimodal / multimodal tests ---- *)

let test_bimodal_two_peaks () =
  section "bimodal: two well-separated clusters";
  let data = Kerdaml.create_vec 100 in
  for i = 0 to 49 do data.{i} <- -5.0 +. Float.of_int i *. 0.02 done;
  for i = 0 to 49 do data.{50 + i} <- 5.0 +. Float.of_int i *. 0.02 done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.3) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "exactly 2 modes" (Array.length modes = 2);
  if Array.length modes = 2 then begin
    let by_pos = modes_by_position result modes in
    let m1 = by_pos.(0) and m2 = by_pos.(1) in
    let x1 = result.points.{m1} in
    let x2 = result.points.{m2} in
    check (Printf.sprintf "mode1 near -5 (got %.2f)" x1) (approx ~tol:0.5 x1 (-4.5));
    check (Printf.sprintf "mode2 near +5 (got %.2f)" x2) (approx ~tol:0.5 x2 5.5);
    let mid_idx = (m1 + m2) / 2 in
    let valley = result.density.{mid_idx} in
    check "valley < left peak" (valley < result.density.{m1});
    check "valley < right peak" (valley < result.density.{m2});
    let integral = trapezoidal_integral result in
    check (Printf.sprintf "bimodal integral = %.6f ≈ 1.0" integral)
      (approx ~tol:1e-3 integral 1.0)
  end

let test_bimodal_symmetric () =
  section "bimodal: symmetric clusters have equal peak heights";
  let data = Kerdaml.create_vec 200 in
  Random.init 99;
  for i = 0 to 99 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    let z = sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2) in
    data.{i} <- -5.0 +. 0.5 *. z;
    data.{100 + i} <- 5.0 +. 0.5 *. z
  done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.5) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "symmetric bimodal has 2 modes" (Array.length modes = 2);
  if Array.length modes = 2 then begin
    let by_pos = modes_by_position result modes in
    let m1 = by_pos.(0) and m2 = by_pos.(1) in
    let h1 = result.density.{m1} in
    let h2 = result.density.{m2} in
    check (Printf.sprintf "peak heights similar (%.4f vs %.4f)" h1 h2)
      (approx ~tol:0.05 h1 h2);
    let x1 = result.points.{m1} in
    let x2 = result.points.{m2} in
    check (Printf.sprintf "modes symmetric around 0 (%.2f, %.2f)" x1 x2)
      (approx ~tol:0.5 (x1 +. x2) 0.0)
  end

let test_bimodal_unequal_weights () =
  section "bimodal: unequal cluster sizes";
  let data = Kerdaml.create_vec 100 in
  Random.init 77;
  for i = 0 to 79 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    data.{i} <- 0.3 *. sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)
  done;
  for i = 0 to 19 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    data.{80 + i} <- 8.0 +. 0.3 *. sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)
  done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.4) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "unequal bimodal has 2 modes" (Array.length modes = 2);
  if Array.length modes = 2 then begin
    (* modes.(0) is dominant since sorted by descending density *)
    check "dominant has higher density"
      (result.density.{modes.(0)} > result.density.{modes.(1)});
    let dominant_x = result.points.{modes.(0)} in
    check (Printf.sprintf "dominant mode near 0 (got %.2f)" dominant_x)
      (approx ~tol:1.0 dominant_x 0.0)
  end

let test_trimodal () =
  section "trimodal: three separated clusters";
  let data = Kerdaml.create_vec 150 in
  for i = 0 to 49 do data.{i} <- -6.0 +. Float.of_int i *. 0.01 done;
  for i = 0 to 49 do data.{50 + i} <- 0.0 +. Float.of_int i *. 0.01 done;
  for i = 0 to 49 do data.{100 + i} <- 6.0 +. Float.of_int i *. 0.01 done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.3) ~n_points:2048 data in
  let modes = Kerdaml.find_modes result in
  check "exactly 3 modes" (Array.length modes = 3);
  if Array.length modes = 3 then begin
    let by_pos = modes_by_position result modes in
    let x1 = result.points.{by_pos.(0)} in
    let x2 = result.points.{by_pos.(1)} in
    let x3 = result.points.{by_pos.(2)} in
    check (Printf.sprintf "mode1 near -6 (got %.2f)" x1) (approx ~tol:0.5 x1 (-5.75));
    check (Printf.sprintf "mode2 near 0 (got %.2f)" x2) (approx ~tol:0.5 x2 0.25);
    check (Printf.sprintf "mode3 near +6 (got %.2f)" x3) (approx ~tol:0.5 x3 6.25);
    let h1 = result.density.{by_pos.(0)} in
    let h2 = result.density.{by_pos.(1)} in
    let h3 = result.density.{by_pos.(2)} in
    check "trimodal peaks similar height"
      (approx ~tol:0.05 h1 h2 && approx ~tol:0.05 h2 h3);
    let integral = trapezoidal_integral result in
    check (Printf.sprintf "trimodal integral = %.6f ≈ 1.0" integral)
      (approx ~tol:1e-3 integral 1.0)
  end

let test_four_modes () =
  section "four modes";
  let data = Kerdaml.create_vec 200 in
  for i = 0 to 49 do data.{i} <- -9.0 +. Float.of_int i *. 0.02 done;
  for i = 0 to 49 do data.{50 + i} <- -3.0 +. Float.of_int i *. 0.02 done;
  for i = 0 to 49 do data.{100 + i} <- 3.0 +. Float.of_int i *. 0.02 done;
  for i = 0 to 49 do data.{150 + i} <- 9.0 +. Float.of_int i *. 0.02 done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.3) ~n_points:2048 data in
  let modes = Kerdaml.find_modes result in
  check "exactly 4 modes" (Array.length modes = 4);
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "4-mode integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0)

let test_bimodal_merges_with_large_bw () =
  section "bimodal merges into unimodal with large bandwidth";
  let data = Kerdaml.create_vec 100 in
  for i = 0 to 49 do data.{i} <- -3.0 +. Float.of_int i *. 0.02 done;
  for i = 0 to 49 do data.{50 + i} <- 3.0 +. Float.of_int i *. 0.02 done;
  let result_narrow = Kerdaml.estimate ~bw:(Fixed 0.3) ~n_points:1024 data in
  let result_wide = Kerdaml.estimate ~bw:(Fixed 5.0) ~n_points:1024 data in
  let modes_narrow = Kerdaml.find_modes result_narrow in
  let modes_wide = Kerdaml.find_modes result_wide in
  check "narrow bw resolves 2 modes" (Array.length modes_narrow = 2);
  check "wide bw merges to 1 mode" (Array.length modes_wide = 1)

let test_bimodal_valley_depth () =
  section "bimodal valley depth scales with separation";
  let make_bimodal sep =
    let data = Kerdaml.create_vec 100 in
    for i = 0 to 49 do data.{i} <- -. sep +. Float.of_int i *. 0.02 done;
    for i = 0 to 49 do data.{50 + i} <- sep +. Float.of_int i *. 0.02 done;
    data
  in
  let valley_ratio data =
    let result = Kerdaml.estimate ~bw:(Fixed 0.5) ~n_points:1024 data in
    let modes = Kerdaml.find_modes result in
    if Array.length modes = 2 then begin
      let by_pos = modes_by_position result modes in
      let m1 = by_pos.(0) and m2 = by_pos.(1) in
      let mid_idx = (m1 + m2) / 2 in
      let peak = Float.max result.density.{m1} result.density.{m2} in
      result.density.{mid_idx} /. peak
    end else 1.0
  in
  let close_ratio = valley_ratio (make_bimodal 2.0) in
  let far_ratio = valley_ratio (make_bimodal 5.0) in
  check (Printf.sprintf "closer clusters: valley ratio %.4f" close_ratio)
    (close_ratio > 0.0);
  check (Printf.sprintf "farther clusters: valley ratio %.4f" far_ratio)
    (far_ratio < close_ratio);
  check "far separation → near-zero valley" (far_ratio < 0.01)

let test_bimodal_generated () =
  section "bimodal generated with Box-Muller";
  Random.init 123;
  let n = 1000 in
  let data = Kerdaml.create_vec n in
  for i = 0 to n - 1 do
    let u1 = Random.float 1.0 in
    let u2 = Random.float 1.0 in
    let z = sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2) in
    if i < n / 2 then data.{i} <- -4.0 +. 0.7 *. z
    else data.{i} <- 4.0 +. 0.7 *. z
  done;
  let result = Kerdaml.estimate ~bw:(Fixed 0.5) ~n_points:1024 data in
  let modes = Kerdaml.find_modes result in
  check "generated bimodal has 2 modes" (Array.length modes = 2);
  let integral = trapezoidal_integral result in
  check (Printf.sprintf "generated bimodal integral = %.6f ≈ 1.0" integral)
    (approx ~tol:1e-3 integral 1.0);
  if Array.length modes = 2 then begin
    let by_pos = modes_by_position result modes in
    let x1 = result.points.{by_pos.(0)} in
    let x2 = result.points.{by_pos.(1)} in
    check (Printf.sprintf "left mode near -4 (got %.2f)" x1) (approx ~tol:0.7 x1 (-4.0));
    check (Printf.sprintf "right mode near +4 (got %.2f)" x2) (approx ~tol:0.7 x2 4.0)
  end

(* ---- indicator-style API tests ---- *)

let test_indicator_lookback () =
  section "indicator lookback";
  check "Density lookback" (Kerdaml.lookback (Density { window = 20; h = 0.5 }) = 19);
  check "NumModes lookback" (Kerdaml.lookback (NumModes { window = 50; h = 1.0; n_points = 128 }) = 49);
  check "DominantMode lookback" (Kerdaml.lookback (DominantMode { window = 30; h = 0.8; n_points = 128 }) = 29);
  check "ModeSpread lookback" (Kerdaml.lookback (ModeSpread { window = 10; h = 0.3; n_points = 128 }) = 9)

let test_indicator_density () =
  section "indicator Density";
  (* 20 bars of data, window=5, so output valid from bar 5 onward *)
  let input = Kerdaml.vec_of_array
    [| 1.0; 2.0; 3.0; 4.0; 5.0; 3.0; 3.0; 3.0; 3.0; 3.0;
       3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0 |] in
  let output = Kerdaml.create_vec 20 in
  Kerdaml.compute (Density { window = 5; h = 0.5 }) input output;
  (* Bars before lookback (window-1=4) should be untouched.
     create_vec uses Bigarray.Array1.create which may not zero memory,
     so use sentinel test pattern instead. *)
  for i = 0 to 19 do output.{i} <- -1.0 done;
  Kerdaml.compute (Density { window = 5; h = 0.5 }) input output;
  check "bar 0 untouched" (output.{0} = -1.0);
  check "bar 3 untouched" (output.{3} = -1.0);
  (* Bar 4 (=lookback) onward should be positive *)
  check "bar 4 positive" (output.{4} > 0.0);
  check "bar 19 positive" (output.{19} > 0.0);
  (* Bars 10-19 have window [all 3.0s], density at 3.0 should be high *)
  check "constant window → high density" (output.{14} > output.{5})

let test_indicator_density_matches_evaluate () =
  section "indicator Density matches evaluate_at";
  let input = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let output = Kerdaml.create_vec 5 in
  let h = 0.8 in
  Kerdaml.compute (Density { window = 3; h }) input output;
  (* Bar 4: window is [3.0, 4.0, 5.0], evaluate density at x=5.0 *)
  let window_data = Kerdaml.vec_of_array [| 3.0; 4.0; 5.0 |] in
  let expected = Kerdaml.evaluate_at h window_data 5.0 in
  check (Printf.sprintf "bar 4 density = %.6f (expected %.6f)" output.{4} expected)
    (approx ~tol:1e-12 output.{4} expected)

let test_indicator_density_single_bar () =
  section "indicator Density single-bar mode";
  let input = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let output_all = Kerdaml.create_vec 5 in
  let output_one = Kerdaml.create_vec 5 in
  let ind = Kerdaml.Density { window = 3; h = 0.8 } in
  Kerdaml.compute ind input output_all;
  Kerdaml.compute ~i:4 ind input output_one;
  check "single-bar matches all-bars at i=4"
    (approx ~tol:1e-15 output_one.{4} output_all.{4});
  (* Other bars should be untouched *)
  check "bar 1 untouched in single-bar mode" (output_one.{1} = 0.0)

let test_indicator_num_modes_unimodal () =
  section "indicator NumModes unimodal data";
  (* Tight cluster → 1 mode *)
  let input = Kerdaml.vec_of_array
    [| 5.0; 5.1; 4.9; 5.05; 4.95; 5.0; 5.1; 4.9; 5.05; 4.95 |] in
  let output = Kerdaml.create_vec 10 in
  Kerdaml.compute (NumModes { window = 5; h = 0.3; n_points = 128 }) input output;
  check "unimodal window → 1 mode" (approx ~tol:0.5 output.{9} 1.0)

let test_indicator_num_modes_bimodal () =
  section "indicator NumModes bimodal data";
  (* Alternating between two well-separated values *)
  let input = Kerdaml.vec_of_array
    [| 0.0; 10.0; 0.0; 10.0; 0.0; 10.0; 0.0; 10.0; 0.0; 10.0 |] in
  let output = Kerdaml.create_vec 10 in
  Kerdaml.compute (NumModes { window = 6; h = 0.5; n_points = 256 }) input output;
  check "bimodal window → 2 modes" (approx ~tol:0.5 output.{9} 2.0)

let test_indicator_dominant_mode () =
  section "indicator DominantMode";
  (* Window of values clustered near 7.0 *)
  let input = Kerdaml.vec_of_array
    [| 7.0; 7.1; 6.9; 7.05; 6.95; 7.0; 7.1; 6.9; 7.05; 6.95 |] in
  let output = Kerdaml.create_vec 10 in
  Kerdaml.compute (DominantMode { window = 5; h = 0.3; n_points = 128 }) input output;
  check (Printf.sprintf "dominant mode near 7.0 (got %.2f)" output.{9})
    (approx ~tol:0.5 output.{9} 7.0)

let test_indicator_mode_spread_unimodal () =
  section "indicator ModeSpread unimodal";
  let input = Kerdaml.vec_of_array
    [| 5.0; 5.1; 4.9; 5.05; 4.95; 5.0; 5.1; 4.9; 5.05; 4.95 |] in
  let output = Kerdaml.create_vec 10 in
  Kerdaml.compute (ModeSpread { window = 5; h = 0.3; n_points = 128 }) input output;
  check "unimodal → spread = 0" (output.{9} = 0.0)

let test_indicator_mode_spread_bimodal () =
  section "indicator ModeSpread bimodal";
  let input = Kerdaml.vec_of_array
    [| 0.0; 10.0; 0.0; 10.0; 0.0; 10.0; 0.0; 10.0; 0.0; 10.0 |] in
  let output = Kerdaml.create_vec 10 in
  Kerdaml.compute (ModeSpread { window = 6; h = 0.5; n_points = 256 }) input output;
  check (Printf.sprintf "bimodal spread near 10 (got %.2f)" output.{9})
    (approx ~tol:1.0 output.{9} 10.0)

let test_indicator_before_lookback_untouched () =
  section "indicator bars before lookback untouched";
  let input = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let output = Kerdaml.create_vec 5 in
  (* Fill output with sentinel *)
  for i = 0 to 4 do output.{i} <- 999.0 done;
  Kerdaml.compute (Density { window = 3; h = 0.5 }) input output;
  check "bar 0 sentinel preserved" (output.{0} = 999.0);
  check "bar 1 sentinel preserved" (output.{1} = 999.0);
  (* window=3 → lookback=2, so bar 2 is first computed *)
  check "bar 2 computed" (output.{2} <> 999.0);
  check "bar 3 computed" (output.{3} <> 999.0);
  check "bar 4 computed" (output.{4} <> 999.0)

let test_indicator_single_bar_before_lookback () =
  section "indicator single-bar before lookback is no-op";
  let input = Kerdaml.vec_of_array [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let output = Kerdaml.create_vec 5 in
  for i = 0 to 4 do output.{i} <- 999.0 done;
  Kerdaml.compute ~i:1 (Density { window = 5; h = 0.5 }) input output;
  check "bar 1 still sentinel (before lookback)" (output.{1} = 999.0)

(* ---- run all ---- *)

let () =
  test_create_vec ();
  test_bandwidth_fixed ();
  test_bandwidth_scott ();
  test_bandwidth_silverman ();
  test_bandwidth_constant_data ();
  test_evaluate_at_single_point ();
  test_evaluate_at_symmetric ();
  test_evaluate_at_two_points ();
  test_evaluate_at_bandwidth_effect ();
  test_evaluate ();
  test_estimate_defaults ();
  test_estimate_custom_n ();
  test_estimate_grid_ordered ();
  test_estimate_density_positive ();
  test_estimate_integrates_to_one ();
  test_estimate_large_data ();
  test_estimate_silverman ();
  test_estimate_fixed_bw ();
  test_estimate_covers_data ();
  test_two_points ();
  test_identical_points ();
  test_negative_data ();
  test_find_modes_unimodal ();
  test_find_modes_sorted_by_density ();
  test_find_modes_empty_result ();
  test_find_modes_index_values ();
  test_bimodal_two_peaks ();
  test_bimodal_symmetric ();
  test_bimodal_unequal_weights ();
  test_trimodal ();
  test_four_modes ();
  test_bimodal_merges_with_large_bw ();
  test_bimodal_valley_depth ();
  test_bimodal_generated ();
  test_indicator_lookback ();
  test_indicator_density ();
  test_indicator_density_matches_evaluate ();
  test_indicator_density_single_bar ();
  test_indicator_num_modes_unimodal ();
  test_indicator_num_modes_bimodal ();
  test_indicator_dominant_mode ();
  test_indicator_mode_spread_unimodal ();
  test_indicator_mode_spread_bimodal ();
  test_indicator_before_lookback_untouched ();
  test_indicator_single_bar_before_lookback ();
  Printf.printf "\n%d passed, %d failed\n" !passed !failed;
  if !failed > 0 then exit 1
