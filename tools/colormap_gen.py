#!/usr/bin/python
'''

Program to generate a Cello-compatible colormap.

Can use either space-delimited RGB values
(such as the ones generated here: http://jdherman.github.io/colormap)

		or

can use matplotlib's built-in colormaps
(see here: http://matplotlib.org/examples/color/colormaps_reference.html)


Note: Might require python 2.x

'''

from matplotlib import cm
from __future__ import print_statement

def get_textfile_colors(filename):
	
	r = []
	g = []
	b = []
	num_lines = 0

	file = open(filename , 'r')
	useful = file.readlines()

	for line in useful:

		red, green, blue  = line.split()
		r.append(float(red)/255.)
		g.append(float(green)/255.)
		b.append(float(blue)/255.)
		num_lines += 1

	file.close()
	
	return r, g, b, num_lines



def get_native_mpl_colors(mapname):

	r = []
	g = []
	b = []
	x = 0
	num_lines = 0

	cmap = cm.get_cmap(mapname)

	for x in range(256): # max colors expected

		r.append(cmap(x)[0])
		g.append(cmap(x)[1])
		b.append(cmap(x)[2])
		num_lines = x

		if(cmap(x) == cmap(x+1)):
			break

	return r, g, b, num_lines



def write_to_file(outfile, reds, greens, blues, length):

	f = open(outfile,'w')

	spaces = "\t\t\t"
	firstline = "# File:\t" + outfile

	f.write(firstline + '\n')
	f.write('\n' + "colormap = [" + '\n')

	for l in range(length-1):
		f.write(spaces + "  ")
		f.write(str(reds[l]) + ", " + str(greens[l]) + \
			    ", " + str(blues[l]) + ',\n')

	f.write(spaces + "  ")
	f.write(str(reds[length-1]) + ", " + str(greens[length-1]) + \
		    ", " + str(blues[length-1]) + '\n')
	f.write(spaces + "];")

	f.close()



promptA = "Would you like to generate a colormap from a " + \
	      "text file (1) or a matplotlib native colormap (2)?"

promptB = "Enter the name of the text file: "

promptC = "Enter the name of the colormap: "

cmType = int(input(promptA))

while ((cmType != 1) and (cmType != 2) and (cmType != 0)):
	print("Invalid argument (0 to exit). ")
	cmType = int(input(promptA))
	
if cmType == 0:

	raise SystemExit

elif cmType == 1:

	infilename = str(raw_input(promptB))

	# maybe add ../input/ to store cmaps in proper folder
	outfilename = "colormap_"

	if infilename.endswith(".txt"):
		outfilename = "colormap_" + infilename[:-4] + ".incl"

	else:
		outfilename = "colormap_" + infilename + ".incl"

	R, G, B, N = get_textfile_colors(infilename)

	write_to_file(outfilename, R, G, B, N)

elif cmType == 2:

	mapname = str(raw_input(promptC))

	outfilename = "colormap_" + mapname + ".incl"

	try:
		R, G, B, N = get_native_mpl_colors(mapname)

	except:
		print("** Invalid. Terminating program. **")
		raise SystemExit

	write_to_file(outfilename, R, G, B, N)
