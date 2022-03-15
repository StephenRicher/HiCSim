find . -name 'simulation.custom.gz' -print0 | 
  while IFS= read -r -d '' file; do 
    count=$(zcat $file | wc -l)
    if [ "${count}" -ne 198789 ]; then
      rm "${file}"
    fi
  done
