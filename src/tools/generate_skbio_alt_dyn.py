#!/usr/bin/env python3

#
# Parse skbio_alt_dyn_impl.hpp
# and generate concrete implementations
# for all inline functions ending in _T
#

#
# Arguments:
#  $1 - namespace
#  $2 - method
#

from skbb_generate_helper import print_header,print_body
import sys

variant = sys.argv[1]
method =sys.argv[2]

#
# ==========================
#

with open('skbio_alt_dyn_impl.hpp','r') as fd:
    lines=fd.readlines()

#
# ==========================
#

# print out the header
print('// Generated from skbio_alt_dyn_impl.hpp (using method %s)'%method);
print('// Do not edit by hand');
print('');

if method in ('direct','indirect',):
    # we are generating skbio_alt_dyn.cpp
    print('#include "skbio_alt_dyn.hpp"');

if method in ('indirect','api',):
    # we referencing the api
    print('#include "skbio_alt_dyn_%s.h"'%variant);

if method in ('direct','api',):
    # we are generating the actual code
    print('#include "skbio_alt_dyn_impl.hpp"');

nmspace = print_header(variant,method)
print_body(method, lines, nmspace)

