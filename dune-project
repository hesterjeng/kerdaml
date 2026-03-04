(lang dune 3.0)

(generate_opam_files true)

(name kerdaml)

(authors "John Hester")

(maintainers "John Hester")

(source
 (github hesterjeng/kerdaml))

(package
 (name kerdaml)
 (synopsis "Univariate kernel density estimation for OCaml")
 (license "GPL-3.0-or-later")
 (homepage "https://github.com/hesterjeng/kerdaml")
 (description "Standalone univariate kernel density estimation library for OCaml")
 (depends ocaml ppx_deriving ppx_hash ppx_compare))
