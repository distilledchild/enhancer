mkdir option1


find . -type f -name '*.png' -print0 |
  while IFS= read -r -d '' file
    do convert "${file}" -resize 1024x "${file%.*}".jpg ;
  done

mv *.jpg ./option1/
