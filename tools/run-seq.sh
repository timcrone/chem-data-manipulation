./gpusimilarity /mnt/fsim/$1
cat mols-to-check | awk -v fname=`basename -s .smi.gz.fsim $1` '{print("python3 python/gpusim_search.py "fname" \047\047 \047"$2"\047 >> "$1".sim")}'
pkill gpusimilarity
