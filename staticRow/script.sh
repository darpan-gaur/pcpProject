# get the name of .cpp files in current directory
for file in *.cpp
do
    # remove .cpp from file name
    name=${file%.cpp}
    echo $name
    # compile the file
    # g++ $file -o $name
    # # run the file

    echo "" > $name.txt
    for i in {1..5}
    do
        n=$((100*i))
        # print n and n/2 in input_par.in
        if [[ $name == *"Batch"* ]]; then
            inpfile="input_batch.in"
        else
            inpfile="input_par.in"
        fi
        echo $n $((n/2))> $inpfile
        echo 0.0 1.0 >> $inpfile
        echo 0.0 0.5 >> $inpfile
        # if name contains Batch echo 16 50 else 16
        if [[ $name == *"Batch"* ]]; then
            echo 16 50 >> $inpfile
        else
            echo 16 >> $inpfile
        fi

        ./$name
        echo -n "$n " >> $name.txt
        for j in {1..5}
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