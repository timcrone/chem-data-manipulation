find -name m\*.tsv -exec awk '{x+=$0;y+=$0^2}END{print "{}", x/NR, sqrt(y/NR-(x/NR)^2)}' {} \;
