# User visible changes in OptimPack


## Version 3.2.0

- Status returned by `BOBYQA`, `COBYLA`, and `NEWUOA` algorithms are now
  enumerations.

- Names of functions, types, and macros no longer start with an underscore
  (reserved to future development of the C language).  Names of types,
  structures, and enumerations no longer have a `_t` suffix.  Names of opaque
  structures have a trailing underscore.  File
  [`tools/update_3p1_to_3p2`](tools/update_3p1_to_3p2) provides a simple Perl
  script to make the necessary changes in your code.

## Version 3.1.0
