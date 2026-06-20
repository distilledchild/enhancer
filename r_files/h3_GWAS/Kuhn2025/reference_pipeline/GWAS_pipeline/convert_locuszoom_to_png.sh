# USAGE: bash convert_locuszoom_to_png.sh <dir with pdfs> <dir for pngs>
# e.g.: bash convert_locuszoom_to_png.sh ./locuszoom_pdfs ./locuszoom_pngs

FILES=$(find $1 -type f -name '*.pdf')
python3 convert_locuszoom_to_png.py $2 $FILES
