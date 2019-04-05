#!/bin/bash

# m(x) excess grid
km_grid="-0.7 -0.5 0.0 1.0 2.0"
count_k_m=5

# m(x) dispersion grid
origin_d_m=0.05
count_d_m=20
last_d_m=0.5

# w(x) excess grid
kw_grid="-0.7 -0.5 0.0 1.0 2.0"
count_k_w=5

# w(x) dispersion grid
origin_d_w=0.05
count_d_w=20
last_d_w=0.5

#accurancy
eps=0.001

m_args="m.txt"
w_args="w.txt"
args="args.txt"
rm -f $m_args $w_args

# WARNING!!! second part was removed due to similar grids for m and w.
# If you want to compute not symmetrical case, please add the second part
echo "Creating arguments..."
ind=0
for km in $km_grid
do
    temp_files[$ind]=`mktemp`
    ./excess $km 1 $km $origin_d_m $count_d_m $last_d_m \
             100001 $eps k ${temp_files[$ind]} &
    pids[$ind]=$!
    ((ind += 1))
done

for pid in ${pids[*]}
do
    wait $pid
done

for temp_file in ${temp_files[*]}
do
    cat $temp_file >> $m_args
    rm -f $temp_file
done

# instead of the second part
cp $m_args $w_args

echo "Collecting them together..."

# $m_args and $w_args files have next form:
#
#   [ km dm s0m s1m ] and [ kw dw s0w s1w ]
#
# We want to get result file in form:
#
#   [ km kw dm dw s0m s1m s0w s1w ]
#
# Note that result file line number is a multiplying of input args file
# line numbers.

rm -f $args

shft_mb=1
shft_me=$count_d_m

for i in `seq $count_k_m`
do
    shft_wb=1
    shft_we=$count_d_w

    for j in `seq $count_k_w`
    do
        sed -n "${shft_mb},${shft_me}p" $m_args | while read line_m
        do
            pars_m=($line_m)
            sed -n "${shft_wb},${shft_we}p" $w_args | while read line_w
            do
                pars_w=($line_w)
                echo "${pars_m[0]} ${pars_w[0]} ${pars_m[1]} ${pars_w[1]}" `
                    `"${pars_m[2]} ${pars_m[3]} ${pars_w[2]} ${pars_w[3]}"`
                    `>> $args
            done
        done

        ((shft_wb += $count_d_w))
        ((shft_we += $count_d_w))
    done

    ((shft_mb += $count_d_m))
    ((shft_me += $count_d_m))
done

rm -f $m_args $w_args
