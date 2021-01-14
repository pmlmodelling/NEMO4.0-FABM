#!/bin/bash

#./inc/clean.sh
#./inc/build.sh

sed -i -e 's#utf8#latin1#'                                 \
       -e 's#\[outputdir=../build\]{minted}#\[\]{minted}#' \
       -e '/graphicspath/ s#{../#{../../#g'                \
          global/packages.tex

cd ./NEMO/main
sed -i -e 's#\\documentclass#%\\documentclass#' -e '/{document}/ s#^#%#'   ../subfiles/*.tex
sed -i    's#\\subfile{#\\input{#'                                         chapters.tex appendices.tex

#latex2html                                   -noimages -local_icons -no_footnode -split 4 -link 2 -dir ../html_LaTeX2HTML $*	NEMO_manual
latex2html -debug -noreuse -init_file ../../l2hconf.pm -local_icons                               -dir ../build/html               NEMO_manual

sed -i -e 's#%\\documentclass#\\documentclass#' -e '/{document}/ s#^%##'   ../subfiles/*.tex
sed -i    's#\\input{#\\subfile{#'                                         chapters.tex appendices.tex
cd -

sed -i -e 's#latin1#utf8#'                                 \
       -e 's#\[\]{minted}#\[outputdir=../build\]{minted}#' \
       -e '/graphicspath/ s#{../../#{../#g'                \
	  global/packages.tex

exit 0
