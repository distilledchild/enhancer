# USAGE: python3 convert_locuszoom_to_png.py <dir for pngs> <any number of pdf paths>
# e.g.: python3 convert_locuszoom_to_png.py ./locuszoom_pngs test1.pdf test2.pdf
# REQUIRES poppler and pdf2image installed (use "conda install poppler")

import sys
from pdf2image import convert_from_path

if len(sys.argv) < 3:
	print("Must pass output directory as first argument, and at least one PDF path")

for path in sys.argv[2:]:
	convert_from_path(path, first_page=1, fmt='png', single_file=True, output_folder=sys.argv[1], output_file=path[:-4])
