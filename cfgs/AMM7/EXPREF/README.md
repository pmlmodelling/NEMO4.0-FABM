# AMM7 Physics Configuration

All inputs for 2005 can currently be found: /work/n01/n01/dapa/NEMO4-TESTS/AMM7/INPUTS

Code can be compiled from the nemo root directory on ARCHER2 using the following:
> module -s restore /work/n01/shared/acc/n01_modules/ucx_env
> CFG=AMM7_TEST
> ARCH=X86_ARCHER2-Cray
> REF=AMM7
> printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -n $CFG -r $REF -m $ARCH -j 0
> ./makenemo -n $CFG -r $REF -m $ARCH -j 16 clean
> ./makenemo -n $CFG -r $REF -m $ARCH -j 16


