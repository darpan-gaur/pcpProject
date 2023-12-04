# get the name of .cpp files in current directory
for file in *.cpp
do
    # remove .cpp from file name
    name=${file%.cpp}
    echo $name
    # compile the file
    g++ $file -o $name
    # # run the file

    echo "" > $name.txt
    # for i in {2, 4, 8, 16, 20, 25, 32}
    for i in 2 4 8 16 20 25 32
    do
        echo $i
        # n=$((100*i))
        # print n and n/2 in input_par.in
        # if name contain Bath then input_batch.in else if name contains seq then input_seq.in else input_par.in
        echo 300 150 > $inpfile
        echo 0.0 1.0 >> $inpfile
        echo 0.0 0.5 >> $inpfile
        # if name contains Batch echo 16 50 else 16
        if [[ $name == *"Batch"* ]]; then
            echo $i 50 >> $inpfile
        else
            echo $i >> $inpfile
        fi

        # ./$name
        # echo -n "$i " >> $name.txt
        for j in {1..1}
        do
            echo "$name $n " 
            ./$name >> $name.txt
        done
    done

    # remove the executable file
    rm $name
done

# g++ oneRow.cpp -o oneRow
# echo "nx, time" > out.txt
# for i in {1..5}
# do
#     n=$((100*i))
#     # print n and n/2 in input_par.in
#     echo $n $((n/2))> input_par.in
#     echo 0.0 1.0 >> input_par.in
#     echo 0.0 0.5 >> input_par.in
#     echo 16 >> input_par.in

#     ./oneRow
#     for j in {1..5}
#     do
#         echo "oneRow $n " 
#         # echo -n "$n " >> out.txt
#         ./oneRow >> out.txt
#     done
#     # echo "" >> out.txt
# done