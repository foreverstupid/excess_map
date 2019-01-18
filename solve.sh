#!/bin/bash

echo "Solving..."

dim=$1
if [ -z $dim ]
then
    dim=1
fi

static_params="-a 0.4 -d 0 -b 1 -s 0.4 -e 7 -p n -r 2 -D $dim -n 1000 -i 300"
plot_data="surface${dim}d.plt"
prog="./neuman"

args=$2
if [ -z $args ]
then
    args="args.txt"
fi

rm -f $plot_data
cat $args | while read line
do
    if [ -z "$line" ]
    then
        echo "" >> $plot_data
    else
        pars=($line)
        echo "$prog $static_params -kK ${pars[4]} ${pars[5]} "`
                                     `"${pars[6]} ${pars[7]}"

        N=`$prog $static_params -kK ${pars[4]} ${pars[5]} \
                                    ${pars[6]} ${pars[7]}`

        if [ $N = "nan" -o $N = "-nan" ]
        then
            echo "NAN"
            N=1.0
        fi

        N=`echo $N | sed -e 's/^[[:space:]]*//'`

        echo "${pars[0]} ${pars[1]} ${pars[2]} ${pars[3]} $N" >> $plot_data
        echo "---> $N"
    fi
done
