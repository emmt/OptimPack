# Changes in OptimPack

This describe the most important (user visible) changes in **OptimPack**.

## Version 3.1.0 (31/10/2019)

This is version 3.1.0 of OptimPack library released on October 31, 2019.

Download: [optimpack-3.1.0.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.1.0/optimpack-3.1.0.tar.bz2)

In this version, the main OptimPack library (*e.g.*, `libopk.so`) no longer
includes Powell's methods which are available in separate libraries as before.
The four libraries provided by OptimPack are **independent**: they can be
linked together in a common executable of plugin without any conflicts, their
functions use specific prefixes.  The following table summarizes their usage.

| Headers          | Link flags     | Prefix    | Description                          |
|:-----------------|:---------------|:----------|:-------------------------------------|
| `<optimpack*.h>` | `-lopk -lm`    | `opk_`    | OptimPack methods (VMLMB, NLCG, ...) |
| `<bobyqa.h>`     | `-lbobyqa -lm` | `bobyqa_` | Powell's BOBYQA method               |
| `<cobyla.h>`     | `-lcobyla -lm` | `cobyla_` | Powell's COBYLA method               |
| `<newuoa.h>`     | `-lnewuoa -lm` | `newuoa_` | Powell's NEWUOA method               |

This version also includes minor changes to fix warnings about unused constants
or uninitialized variables.  The code can be compiled with flags `-Wall
-Werror`.

## Version 3.0.1 (20/04/2017)

This version has been released to work with the new Julia interface of
OptimPack.  Few things have changed:
- Build separate libraries.
- Code cleanup in NEWUOA.


## Version 3.0.0

This is version 3.0.0 of OptimPack library released on the 15th of March 2017.

Download: [optimpack-3.0.0.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.0.0b/optimpack-3.0.0.tar.gz)

Most important changes:

- Major update of the **VMLM-B** optimizer for bound constraints optimization.
  This optimizer can be used without bounds and is then similar to Nocedal's
  **VMLM** or Liu & Nocedal's **L-BFGS** algorithms.  It can also emulates the
  Benson & Moré's **BLMVM** algorithm.  This optimizer supersedes (and
  replaces) all quasi-Newton (variable metric) algorithms that were implemeted
  in OptimPack.


## Version 3.0.0b

This is beta pre-3.0.0 version of OptimPack library released on the 19th of
December 2015.

Download: [optimpack-3.0.0b.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.0.0b/optimpack-3.0.0b.tar.gz)

Most important changes:

- Line search can be chosen to run CUTEst tests.

- General VMLMN optimizer to emulate VMLM, L-BFGS, BLMVM and VMLM-B. This
  optimizer will be renamed as VMLMB and will replace all other quasi-Newton
  (variable metric) algorithms.

- The bound constraints are simply specified when the optimizer is created.


## Version 3.0.0a

This is alpha pre-3.0.0 version of OptimPack library released on the 19th of
December 2015.

Download: [optimpack-3.0.0a.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.0.0b/optimpack-3.0.0a.tar.gz)

Most important changes:

- API has changed (the various optimizers can be used in a more uniform way).

- New optimizer VMLMB for bound constrained optimization.

- CUTEst interface to test the algorithms on various problems.


## Version 2.0.1

This is version 2.0.1 of OptimPack released on 10 Mar 2015.

Download: [optimpack-2.0.1.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.0.0b/optimpack-2.0.1.tar.gz)


## Version 2.0.0

This is version 2.0.0 of OptimPack released on 10 Mar 2015.

Download: [optimpack-2.0.0.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.0.0b/optimpack-2.0.0.tar.gz)
