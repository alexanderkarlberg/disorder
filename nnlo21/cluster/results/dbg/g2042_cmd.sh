#!/bin/bash
export P2BSLICE=1 P2BDEBUG=1 VEGAS_EQUAL=2
exec /ptmp/mpp/akarlber/disorder-nnlo21/bin3/nlo31 r 1000000 6 2042 1d-9 2 0 0 psmc 0 0 1d-12
