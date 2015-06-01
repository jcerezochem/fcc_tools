#!/bin/bash

i=0

cd tests/
echo "****************************"
(( i++ ))
echo "TEST $i: GAMESS"
echo "****************************"
echo "gen_fcc_state -i HQ_gamess.out -ft gms"
../src/state_generator/gen_fcc_state -i HQ_gamess.out -ft gms
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (LOG)"
echo "****************************"
echo "gen_fcc_state -i HQ_gaussian.log"
../src/state_generator/gen_fcc_state -i HQ_gaussian.log
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: PSI4"
echo "****************************"
echo "gen_fcc_state -i h2o_psi4.out -ft psi4"
../src/state_generator/gen_fcc_state -i h2o_psi4.out -ft psi4
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: MOLCAS"
echo "****************************"
echo "gen_fcc_state -i s0_mq-GP_edit.UnSym -ft molcas"
../src/state_generator/gen_fcc_state -i s0_mq-GP_edit.UnSym -ft molcas
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: MOLPRO"
echo "****************************"
echo "gen_fcc_state -i h2o_molpro.out -ft molpro"
../src/state_generator/gen_fcc_state -i h2o_molpro.out -ft molpro
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (FCHK) -Linear molecule"
echo "****************************"
echo "gen_fcc_state -i diacetylene.fchk"
../src/state_generator/gen_fcc_state -i diacetylene.fchk
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (LOG) -Linear molecule"
echo "****************************"
echo "gen_fcc_state -i diacetylene.log"
../src/state_generator/gen_fcc_state -i diacetylene.log
echo ""

cd ..

