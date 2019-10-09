export name=RNAP_XLs

cp ../GSMs_cl-1/sample_A.txt Scores_A.txt
cp ../GSMs_cl-1/sample_B.txt Scores_B.txt

python ../scripts/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path ../GSM_cl-1 --mode cuda --align --density density.txt --gridsize 3.0 --gnuplot --scoreA Scores_A.txt --scoreB Scores_B.txt > clustering.log &
