#!/bin/bash

i=0

cd tests/
echo "++++++++++++++++++++++++++++++"
echo " ++++++++++++++++++++++++++++ "
echo "       gen_fcc_state          "
echo " ++++++++++++++++++++++++++++ "
echo "++++++++++++++++++++++++++++++"
echo " "
echo "****************************"
(( i++ ))
echo "TEST $i: GAMESS"
echo "****************************"
echo "gen_fcc_state -i HQ_gamess.out -ft gms"
../generators/gen_fcc_state -i HQ_gamess.out -ft gms
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (LOG)"
echo "****************************"
echo "gen_fcc_state -i HQ_gaussian.log"
../generators/gen_fcc_state -i HQ_gaussian.log
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: ORCA (HESS)"
echo "****************************"
echo "gen_fcc_state -ft orca -i HQ_orca.out"
../generators/gen_fcc_state -ft orca -i HQ_orca.hess
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: PSI4"
echo "****************************"
echo "gen_fcc_state -i h2o_psi4.out -ft psi4"
../generators/gen_fcc_state -i h2o_psi4.out -ft psi4
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: MOLCAS"
echo "****************************"
echo "gen_fcc_state -i s0_mq-GP_edit.UnSym -ft molcas"
../generators/gen_fcc_state -i s0_mq-GP_edit.UnSym -ft molcas
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: MOLPRO"
echo "****************************"
echo "gen_fcc_state -i h2o_molpro.out -ft molpro"
../generators/gen_fcc_state -i h2o_molpro.out -ft molpro
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (FCHK) -Linear molecule"
echo "****************************"
echo "gen_fcc_state -i diacetylene.fchk"
../generators/gen_fcc_state -i diacetylene.fchk
echo ""

echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (LOG) -Linear molecule"
echo "****************************"
echo "gen_fcc_state -i diacetylene.log"
../generators/gen_fcc_state -i diacetylene.log
echo ""
echo ""

echo "++++++++++++++++++++++++++++++"
echo " ++++++++++++++++++++++++++++ "
echo "      gen_fcc_dipfile         "
echo " ++++++++++++++++++++++++++++ "
echo "++++++++++++++++++++++++++++++"
echo " "
echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (FCHK) "
echo "****************************"
echo "gen_fcc_dipfile -i azabenz-ccl4-s2-neq-fr.fchk"
../generators/gen_fcc_dipfile -i azabenz-ccl4-s2-neq-fr.fchk
echo ""
echo "****************************"
(( i++ ))
echo "TEST $i: GAUSSIAN (FCHK) "
echo "****************************"
echo "gen_fcc_dipfile -i azabenz-ccl4-s2-neq-fr.log"
../generators/gen_fcc_dipfile -i azabenz-ccl4-s2-neq-fr.log
echo ""
echo "****************************"
(( i++ ))
echo "TEST $i: PSI4 "
echo "****************************"
echo "gen_fcc_dipfile -i s-methyloxirane_gs_aug.out -ft psi4 -Si 1 -Sf 2 -noder"
../generators/gen_fcc_dipfile -i s-methyloxirane_gs_aug.out -ft psi4 -Si 1 -Sf 2 -noder
echo ""

cd ..

