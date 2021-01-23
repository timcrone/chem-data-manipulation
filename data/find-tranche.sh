curl http://zinc15.docking.org/substances/$1/ -o /tmp/x > /dev/null 2>&1 ; grep  tranche_name /tmp/x | awk -F= '{print $3}' | awk -F'"' '{print $1}'
