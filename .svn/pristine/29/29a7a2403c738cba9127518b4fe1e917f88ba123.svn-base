#!/bin/bash

./inc/clean.sh
./inc/build.sh

mkdir -p html_htlatex
cd tex_main
htlatex NEMO_manual "NEMO_manual,2" "" "-d../html_htlatex/" "-shell-escape"
cd -

exit 0
