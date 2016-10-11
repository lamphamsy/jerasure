#!/bin/bash -e

trap "rm -fr T Coding"  EXIT
KB=1024
MB=$KB*$KB

k=12
m=6
w=8
n=1000
pk_size=$MB

s1=10*$MB
s2=2
dd if=/dev/urandom of=T bs=$s1 count=$s2

echo "*-----------PARAMETERS--------------*"
echo "Erasure codes params: k = $k, m = $m"
echo "Galois field over GF(2^$w)"
echo "Object size: $(($s1*$s2/($MB))) (MB)"
echo "Number of encodings: $n"

echo "*-----------ENCODING USING GF MATRIX----------------*"
echo "Reed-Solomon"
./encoder T $k $m reed_sol_van $w 0 0 $n

echo "reed_sol_r6_op"
./encoder T $k $m reed_sol_van $w 0 0 $n

echo "OUR CODE"
./encoder T $k $m ec_pcm $w 0 0 $n

echo "*-----------ENCODING USING BINARY MATRIX----------------*"
echo "cauchy_orig"
./encoder T $k $m cauchy_orig $w $pk_size  0 $n

echo "cauchy_good"
./encoder T $k $m cauchy_good $w $pk_size  0 $n

echo "OUR CODE"
./encoder T $k $m ec_pcm $w $pk_size  0 $n cauchy
