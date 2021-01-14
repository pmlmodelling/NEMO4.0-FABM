#!/bin/sh

clean() {
    printf "\t¤ Clean previous build"
    find latex/$1/build -mindepth 1 -delete

    echo
}

build() {
    printf "\t¤ Generation of the PDF format\n"
    latexmk -r  ./latex/global/latexmk.pl -pdfxe ./latex/$1/main/$1_manual \
#	1> /dev/null
    [ -f ./latex/$1/build/$1_manual.pdf ] && mv ./latex/$1/build/$1_manual.pdf .
    echo
}
