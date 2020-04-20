cp ../data/covid-clean-04192046.sim.gz* ../data/covid-clean-04192046.sim.gz ; zcat ../data/covid-clean-04192046.sim.gz | awk '{delta = $4 > $3 ? $4 - $3 : $3 - $4 ; if (delta > 0.4)  print $0}'
