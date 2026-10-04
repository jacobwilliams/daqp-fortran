---
project: daqp-fortran
summary: Modern Fortran port of DAQP, a dual active-set solver for dense convex quadratic programs
author: Jacob Williams
github: https://github.com/jacobwilliams
project_github: https://github.com/jacobwilliams/daqp-fortran
license: MIT
src_dir: ./src
output_dir: ./doc
preprocess: true
preprocessor: gfortran -E -cpp
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
fpp_extensions: F90
display: public
         protected
         private
source: true
graph: true
search: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            ieee_arithmetic:https://gcc.gnu.org/onlinedocs/gfortran/IEEE-modules.html
---

{!README.md!}
