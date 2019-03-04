#!/bin/bash

echo "Initialization..."

# variables settings

dim=$1
if [ -z $dim ]
then
    dim=1
fi

static_params="-a 0.4 -d 0 -b 1 -s 0.4 -e 7 -p n -r n -D $dim -n 10000 -i 300"
plot_data="surface${dim}d.plt"
args="args.txt"

threads=$2
if [ -z $threads ]
then
    threads=4
fi
(( part_size=`wc $args | awk '{print $1}'` / $threads ))
(( threads -= 1 ))
thread_suf_len=`expr length $threads`
temp_dir=`mktemp -d`
prefix="$temp_dir/args"
prog="./neuman"



# main claculating logic

split $args -d -a $thread_suf_len $prefix -l $part_size 

echo "Solving..."

pid_counter=0
for i in `seq -w 0 $threads`
do
    args="$prefix$i"
    plot_part="$temp_dir/$i"
    rm -f $plot_part

    cat $args | while read line
    do
        if [ -z "$line" ]
        then
            echo "" >> $plot_part
        else
            pars=($line)
            N=`$prog $static_params -kK ${pars[4]} ${pars[5]} \
                                        ${pars[6]} ${pars[7]}`

            if [ $N = "nan" -o $N = "-nan" ]
            then
                echo "NAN"
                N=1.0
            fi

            N=`echo $N | sed -e 's/^[[:space:]]*//'`

            echo "${pars[0]} ${pars[1]} ${pars[2]} ${pars[3]} $N" `
               ` >> $plot_part
            echo "Solved: $N (${pars[4]} ${pars[5]} ${pars[6]} ${pars[7]})"
        fi
    done &

    pids[$pid_counter]=$!
    (( pid_counter += 1 ))
done



# waiting for finishing

for pid in ${pids[*]}
do
    wait $pid
done



# colleting parts together

echo
echo "Collecting..."

rm -f $plot_data
for i in `seq -w 0 $threads`
do
    args="$prefix$i"
    plot_part="$temp_dir/$i"

    cat $plot_part >> $plot_data
done
rm -rf $temp_dir

echo "Done"
