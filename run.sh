#!/bin/bash

export OMP_NUM_THREADS=$4

./numericalIntegration $1 $2 $3
