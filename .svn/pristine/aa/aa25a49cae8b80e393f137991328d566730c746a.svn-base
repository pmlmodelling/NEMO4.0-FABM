#!/bin/sh

#set -evx

if [[ $* != '' ]]; then
	 if [[ $1 = 'all' ]]; then models='NEMO SI3 TOP'; else models=$1; fi
else
	 models='NEMO'
fi

for model in $models; do
	 [[ $model =~ ^(SI3|TOP)$ ]] && continue

    files=$( find latex/$model -type f -name *.tex )

	 echo '¤ Missing namelist groups in '$model' manual'

	 for nlst in $( ls namelists ); do
        [[ ! $( grep \\nlst{$nlst} $files ) ]] && printf '%s ' $nlst
	 done

    echo; echo
	 echo '¤ Vanished index entries in '$model' manual (\{hf,jp,key,mdl,ngn,nlst,np,rou}{...})'

    for file in $files; do

        items=$( grep -Eho "(hf|jp|key|mdl|ngn|nlst|np|rou){[a-zA-Z0-9_\]*}" $file | sort -u )

        if [[ $items == '' ]]; then
            continue
        else
            printf ${file/latex\/*\/subfiles\/}': '
        fi

        for item in $items; do

            arg=$( echo $item | sed 's/.*{\([^}]*\)}/\1/' | tr -d '\\' )

		      if [[ ( $item =~ ^hf   && ! $( find ../src -type f -name   $arg.h90               ) ) || \
                  ( $item =~ ^jp   && ! $( grep ":: *$arg"             ../src/OCE/par_oce.F90 ) ) || \
                  ( $item =~ ^key  && ! $( grep -ri "#if .* $arg"      ../src                 ) ) || \
		            ( $item =~ ^mdl  && ! $( find ../src -type f -name   $arg.[Ff]90            ) ) || \
                  ( $item =~ ^ngn  && ! $( grep \&$arg                 namelists/*            ) ) || \
                  ( $item =~ ^nlst && ! -f namelists/$arg                                       ) || \
                  ( $item =~ ^np   && ! $( grep " $arg *="             namelists/*            ) ) || \
                  ( $item =~ ^rou  && ! $( grep -ri "SUBROUTINE *$arg" ../src                 ) )      ]]; then
                printf $item' '
            fi

        done

        echo

    done

done

echo
echo '¤ Namelist parameters unfollowing naming conventions (^[cdlnr]n_* or uppercase somewhere)'

for nlst in $( ls namelists ); do
    np_list=$( sed '/^ *[!/&]/d; s|[!/&].*||' namelists/$nlst | tr -d ' ' | cut -d= -f1 )
    array=()

    for np in ${np_list}; do

        if [[ ! ${np:0:3} =~ ^[cdlnr]n_$ || $( echo $np | grep [A-Z] ) ]]; then
            array+=$np' '
        fi

    done

	 if (( ${#array[@]} != 0 )); then
        printf '%15s: ' $nlst
        echo ${array[@]}
    fi

done

exit 0
