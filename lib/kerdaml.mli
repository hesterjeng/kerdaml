type vec = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

type bandwidth =
  | Scott
  | Silverman
  | Fixed of float

type estimate = {
  points : vec;
  density : vec;
}

val create_vec : int -> vec
val vec_of_array : float array -> vec

val compute_bandwidth : bandwidth -> vec -> float
(** Compute the bandwidth value for the given data. *)

val evaluate_at : float -> vec -> float -> float
(** [evaluate_at h data x] — evaluate KDE with bandwidth [h] at point [x]. *)

val estimate : ?bw:bandwidth -> ?n_points:int -> vec -> estimate
(** Estimate density over a uniform grid. Defaults: Scott bandwidth, 512 points.
    Returns a record with [points] (grid) and [density] (estimates). *)

val evaluate : ?bw:bandwidth -> vec -> float -> float
(** Convenience: compute bandwidth then evaluate at a single point. *)

val find_modes : estimate -> int array
(** Return indices of local maxima in the density curve, sorted by
    descending density. Use [result.points.{i}] and [result.density.{i}]
    to get the x-value and density at each mode. *)

(** {1 Indicator-style API}

    Per-bar scalar outputs computed from a lookback window of prices.
    Designed to match tacaml's pattern:
    - [?i] for single-bar computation (live trading)
    - Input/output are vecs (bigarray rows from a data matrix)
    - Writes in place, no allocation visible to caller *)

type indicator =
  | Density of { window : int; h : float }
  | NumModes of { window : int; h : float; n_points : int }
  | DominantMode of { window : int; h : float; n_points : int }
  | ModeSpread of { window : int; h : float; n_points : int }
[@@deriving show, eq, hash, compare]
(** KDE-derived indicators. [window] is the lookback size. [h] is bandwidth.
    - [Density]: density at current price (how "typical" it is)
    - [NumModes]: number of modes in the window (1 = unimodal, 2+ = multimodal)
    - [DominantMode]: x-value of the highest-density mode
    - [ModeSpread]: distance between lowest and highest mode (0 if unimodal) *)

val lookback : indicator -> int
(** Number of past bars required before output is valid. *)

val compute : ?i:int -> indicator -> vec -> vec -> unit
(** [compute indicator input output] — for each bar [i] (or just bar [i] if
    [~i] is given), look back [window] bars in [input], compute the KDE-derived
    scalar, and write it to [output.{i}]. Bars before [lookback] are left
    untouched. *)
