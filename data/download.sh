echo $1 $(curl http://zinc15.docking.org/substances/$1/ 2>/dev/null | grep tranche_name | awk -F= '{print $3}' | awk -F'"' '{print $1}')
