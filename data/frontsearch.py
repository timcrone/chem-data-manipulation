import sys

terms = None
with open("covid.search", "r") as f:
  terms = f.readlines()
with open("/mnt/ssd-data/files.docking.org/3D/alldata.txt") as f:
  l = f.readline()
  count = 0
  found_count = 0
  while l and terms:
    for index, m in enumerate(terms):
      term = m.strip()
      if l.startswith(term):
         print(l)
         terms.pop(index)
         found_count += 1
         print(f"\nfound: {found_count}", file=sys.stderr)
         continue
    count += 1
    print(count, file=sys.stderr, end="\r")
    l = f.readline()
if terms:
  print(f"not found: {terms}")
