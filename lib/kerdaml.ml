type vec = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

type bandwidth =
  | Scott
  | Silverman
  | Fixed of float

type estimate = {
  points : vec;
  density : vec;
}

let create_vec n = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout n

let vec_of_array a =
  let n = Array.length a in
  let v = create_vec n in
  for i = 0 to n - 1 do
    v.{i} <- Array.unsafe_get a i
  done;
  v

let vec_length v = Bigarray.Array1.dim v

(* -- stats helpers -- *)

let mean v =
  let n = vec_length v in
  let s = ref 0.0 in
  for i = 0 to n - 1 do
    s := !s +. v.{i}
  done;
  !s /. Float.of_int n

let std v =
  let n = vec_length v in
  let m = mean v in
  let s = ref 0.0 in
  for i = 0 to n - 1 do
    let d = v.{i} -. m in
    s := !s +. d *. d
  done;
  sqrt (!s /. Float.of_int (n - 1))

let quantile v q =
  (* sort a temporary float array copy *)
  let n = vec_length v in
  let tmp = Array.init n (fun i -> v.{i}) in
  Array.sort Float.compare tmp;
  let idx = q *. Float.of_int (n - 1) in
  let lo = int_of_float (floor idx) in
  let lo = max 0 (min lo (n - 2)) in
  let hi = lo + 1 in
  if lo = hi then tmp.(lo)
  else
    let frac = idx -. Float.of_int lo in
    tmp.(lo) *. (1.0 -. frac) +. tmp.(hi) *. frac

(* -- kernel -- *)

let inv_sqrt_2pi = 1.0 /. sqrt (2.0 *. Float.pi)

let[@inline] gaussian_kernel u = inv_sqrt_2pi *. exp (-0.5 *. u *. u)

(* -- bandwidth -- *)

let compute_bandwidth bw data =
  let n = vec_length data in
  let nf = Float.of_int n in
  match bw with
  | Fixed h -> h
  | Scott ->
    let sigma = std data in
    1.06 *. sigma *. (nf ** (-0.2))
  | Silverman ->
    let sigma = std data in
    let iqr = quantile data 0.75 -. quantile data 0.25 in
    0.9 *. (Float.min sigma (iqr /. 1.34)) *. (nf ** (-0.2))

(* -- evaluation -- *)

let evaluate_at h data x =
  let n = vec_length data in
  let inv_h = 1.0 /. h in
  let s = ref 0.0 in
  for i = 0 to n - 1 do
    let u = (x -. data.{i}) *. inv_h in
    s := !s +. gaussian_kernel u
  done;
  !s *. inv_h /. Float.of_int n

let evaluate ?(bw = Scott) data x =
  let h = compute_bandwidth bw data in
  evaluate_at h data x

let estimate ?(bw = Scott) ?(n_points = 512) data =
  let h = compute_bandwidth bw data in
  (* find min/max *)
  let lo = ref data.{0} in
  let hi = ref data.{0} in
  let n = vec_length data in
  for i = 1 to n - 1 do
    let v = data.{i} in
    if v < !lo then lo := v;
    if v > !hi then hi := v
  done;
  let pad = 3.0 *. h in
  let x_lo = !lo -. pad in
  let x_hi = !hi +. pad in
  let step = (x_hi -. x_lo) /. Float.of_int (n_points - 1) in
  let points = create_vec n_points in
  let density = create_vec n_points in
  for i = 0 to n_points - 1 do
    let x = x_lo +. Float.of_int i *. step in
    points.{i} <- x;
    density.{i} <- evaluate_at h data x
  done;
  { points; density }

let find_modes result =
  let n = vec_length result.density in
  (* Use >= on the left to handle plateau peaks (symmetric data can produce
     two adjacent grid points with identical density). This reports the
     first index of each plateau as the mode. *)
  let count = ref 0 in
  for i = 1 to n - 2 do
    if result.density.{i} >= result.density.{i - 1}
       && result.density.{i} > result.density.{i + 1}
    then incr count
  done;
  let modes = Array.make !count 0 in
  let k = ref 0 in
  for i = 1 to n - 2 do
    if result.density.{i} >= result.density.{i - 1}
       && result.density.{i} > result.density.{i + 1}
    then (modes.(!k) <- i; incr k)
  done;
  Array.sort (fun a b ->
    Float.compare result.density.{b} result.density.{a}
  ) modes;
  modes

(* -- indicator-style API -- *)

let hash_fold_int = Ppx_hash_lib.Std.Hash.fold_int
let hash_fold_float = Ppx_hash_lib.Std.Hash.fold_float
let compare_int = Int.compare
let compare_float = Float.compare

type indicator =
  | Density of { window : int; h : float }
  | NumModes of { window : int; h : float; n_points : int }
  | DominantMode of { window : int; h : float; n_points : int }
  | ModeSpread of { window : int; h : float; n_points : int }
[@@deriving show, eq, hash, compare]

let lookback ind =
  let window = match ind with
    | Density { window; _ }
    | NumModes { window; _ }
    | DominantMode { window; _ }
    | ModeSpread { window; _ } -> window
  in
  window - 1

let compute_bar_density h window input bar_idx =
  let x = input.{bar_idx} in
  let start = bar_idx - window + 1 in
  let inv_h = 1.0 /. h in
  let s = ref 0.0 in
  for j = start to bar_idx do
    let u = (x -. input.{j}) *. inv_h in
    s := !s +. gaussian_kernel u
  done;
  !s *. inv_h /. Float.of_int window

let compute_bar_with_estimate h n_points window input bar_idx =
  let start = bar_idx - window + 1 in
  let win = create_vec window in
  for j = 0 to window - 1 do
    win.{j} <- input.{start + j}
  done;
  estimate ~bw:(Fixed h) ~n_points win

let compute_bar ind input bar_idx =
  match ind with
  | Density { window; h } ->
    compute_bar_density h window input bar_idx
  | NumModes { window; h; n_points } ->
    let est = compute_bar_with_estimate h n_points window input bar_idx in
    let modes = find_modes est in
    Float.of_int (Array.length modes)
  | DominantMode { window; h; n_points } ->
    let est = compute_bar_with_estimate h n_points window input bar_idx in
    let modes = find_modes est in
    if Array.length modes > 0 then est.points.{modes.(0)}
    else Float.nan
  | ModeSpread { window; h; n_points } ->
    let est = compute_bar_with_estimate h n_points window input bar_idx in
    let modes = find_modes est in
    let n = Array.length modes in
    if n < 2 then 0.0
    else begin
      let lo = ref Float.infinity in
      let hi = ref Float.neg_infinity in
      for k = 0 to n - 1 do
        let x = est.points.{modes.(k)} in
        if x < !lo then lo := x;
        if x > !hi then hi := x
      done;
      !hi -. !lo
    end

let compute ?i ind input output =
  let lb = lookback ind in
  let len = vec_length input in
  match i with
  | Some bar_idx ->
    if bar_idx >= lb then
      output.{bar_idx} <- compute_bar ind input bar_idx
  | None ->
    for bar_idx = lb to len - 1 do
      output.{bar_idx} <- compute_bar ind input bar_idx
    done
