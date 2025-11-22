#!/bin/bash

gcc -std=c11 -O2 ukf.c ukf_test.c -o ukf_test -lm -lc
