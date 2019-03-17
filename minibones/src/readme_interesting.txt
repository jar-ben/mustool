Interesting Variables
---------------------
One can specify a set of variables for which backbones should be computed.
This is specified by a 'b-line' in the input DIMACS. The line starts with the letter 'b' and is
followed by number of variables that should be considered for backbone computation. The line *must*
terminate by '0'.
Currently works for core-based approach (-e) and chunk size 1, lower bound (default) approach.

--Mikolas
20 Feb 2012
