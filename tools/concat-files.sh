find . -maxdepth 3 -mindepth 3 -type f -name $1????.txt -exec tail -n +2 {} \; > $1.tsv
