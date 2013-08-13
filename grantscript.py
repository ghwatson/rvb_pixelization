#!/usr/bin/python

# Derived from Ann's script (annscript.py).

##### don't forget to change the directory
#- build directory
buildDir = "../../build"

#- run directory name
runDir = "square7"

#- maximum number of sites to add in a ratio step
max_add = 8

#- the first simulation size
min_L= 56

#- last (or > last) simulation size
max_L = 73

#- Step size for system sizes (LxL simulations)
step_L = 4

#- length of simulations (a million MC steps times this number)
mult = 40

# Renyi = 2 #S_\alpha
# Lattice = 4 #3=triang, 4=sqr, 6=honeycomb

## -- Remnants of Mario's script --
#parameter["bash"]  = "-shift shift.txt -mult 1" 
#-shift specifies the shift-file (since shift.txt is generated automatically leave it as is)
#parameter["args"]   = [["L", "H"], "g"]
#-mult determines how many samples we want: mult == 1 <--> thermalization = 0.1M, samples = 1M
#parameter["sq"]	= "-r 30h --mpp 1g"
#sharcnet specific commands, see/change lauch_programm
#parameter["files"] = ["../../build/examples/sim"]
#where the executable is located
#parameter["cmake"] = "-DUSE_S:STRING=2 -DUSE_GRID:STRING=4"
#USE_S is the renyi index and USE_GRID is the grid type (3=tri, 4=sqr, 6=hex)
#if you change this you need to recompile
## -- End of Mario's script leftovers --

import math
import os
import subprocess

# Returns a list of tuples along a line between the beginning and end.
# Only works for horizontal lines as of now. The points are ordered from
# right-to-left.
def get_horizontal_line(begin,end):
	line = []
	if begin[0] < end[0]:
		sgn = 1
	else:
		sgn = -1
	for i in xrange( sgn*(end[0] - begin[0]) + 1 ):
		newpt = (begin[0] + i*sgn, begin[1])
		line.append(newpt)
	return line
	
# Provides the rasterized version of a disk using a slightly altered version
# of the Bresenham circle algorithm:
#	 http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
# The ordering of the points in the set is like so:
#			  9
#		   8  7  6
#		5  4  3  2  1
#		  10 11  12
#			 13
def raster_disk(x0,y0,radius):
	# Generate the coordinates.
	f = 1 - radius
	ddF_x = 1
	ddF_y = -2 * radius
	x = 0
	y = radius
	
	coords = []
	
	coords.append((x0, y0 + radius))
	coords.append((x0, y0 - radius))
	right = (x0 + radius, y0)
	left = (x0 - radius, y0)
	coords = coords + get_horizontal_line(right, left)
	
	while x < y :
#		 ddF_x == 2 * x + 1
#		 ddF_y == -2 * y
#		 f == x*x + y*y - radius*radius + 2*x - y + 1
		if f >= 0:
			y -= 1
			ddF_y += 2
			f += ddF_y
		x += 1
		ddF_x += 2
		f += ddF_x	
		
		line1 = get_horizontal_line((x0 + x, y0 + y), (x0 - x, y0 + y))
		line2 = get_horizontal_line((x0 + x, y0 - y), (x0 - x, y0 - y))
		line3 = get_horizontal_line((x0 + y, y0 + x), (x0 - y, y0 + x))
		line4 = get_horizontal_line((x0 + y, y0 - x), (x0 - y, y0 - x))
		
		coords = coords + line1 + line2 + line3 + line4

	return coords

# Prints matrix so that x,y is oriented in the conventional manner.
# Looks like:
# y
# ^
# |
# ---> x
def printRegion(mat):
	for bla in reversed(mat):
		print ' '.join(map(str, bla))
	print('\n')

# Assumes that the elements of the matrix are string-able.
# Maps like:
# ---> y
# |
# v
# x
def printMatToFile(mat,filename,option):
	m = open(filename, option)
	for i in mat:
		i = [str(j) for j in i]
		m.write(" ".join(i) + "\n") 
	m.write('\n')
	m.close()
	
class RegionData(object):
	def __init__(self,L,create_layers,create_bin_partition):
		# Get the layers (ex shells of disk, or columns of square).
		self.layers, self.points = create_layers(L)
		self.values = [0 for x in xrange(self.layers[-1])]
		
		self.bins = create_bin_partition(self.layers)
	
	# Get the rth bin of values.
	def get_bin(self, r):
		if r is not 0:
			begin = self.bins[r-1]
			end = self.bins[r]
		else:
			begin = end = 0
		return self.values[begin:end]
	
	# Set the rth bin to value.
	def set_bin(self, r, value):
		if r is not 0:
			begin = self.bins[r-1]
			end = self.bins[r]
		else:
			begin = 0
			end = 1
		self.values[begin:end] = [value for x in xrange(end-begin)]
		
	# Returns a 2darray of 0s and 1s that represents the current cartesian 
	# state of the data.
	def generate_shiftmap(self,canvas_L):
		canvas = [[0 for x in xrange(canvas_L)] for x in xrange(canvas_L)]
		for i, value in enumerate(self.values):
			# Break once a 0 is reached.
			if not value:
				break
			# Set the corresponding point of value to 1.
			(x,y) = self.points[i]
			canvas[x][y] = value
		return canvas

#removes duplicates in a list whilst preserving order
def remove_duplicates(seq):
	seen = set()
	seen_add = seen.add
	return [ x for x in seq if x not in seen and not seen_add(x) ]

# Create circles with integer radii.
def create_layers(L):
	layer_idx = []
	points = []
	x0 = y0 = int(math.floor(L/4))
	max_r = int(math.ceil(L/4))
	
	# Initialize first layer.
 	points.append((x0,y0))
 	layer_idx.append(1)

	# Collect all the disks full of points.
	rasterdisks = {}
	rasterdisks[0] = [points[0]]
	for r in xrange(1,max_r):
		rasterdisks[r] = raster_disk(x0,y0,r)

	# Build the layers.
	for r in xrange(1,max_r):
   		layer_of_points = [x for x in rasterdisks[r] if x not in rasterdisks[r-1]]
		layer_of_points = remove_duplicates(layer_of_points)
		layer_of_points = order_disk_layer(layer_of_points,r,(x0,y0))
		layer_of_points = remove_duplicates(layer_of_points)
		
		points = points + layer_of_points
		layer_idx.append(len(points))
		
	return layer_idx, points

# Order the points in a layer starting from (R,0) and going around CCW.
# TODO: we only needed the octant here for ordering...consider just working
# with octants as opposed to the full thing.
# Note: The diagonal values get double-counted.
def order_disk_layer(layer,R,origin):
	(x0,y0) = origin
	# Move points to origin
	layer = [(pt[0] - x0,pt[1] - y0) for pt in layer]
	
	# Remove points not in first octant.
	quad0 = [ pt for pt in layer if (pt[0] > 0 and pt[1] > 0) ]
		
	# Perform sorting.
	quad0 = sorted(quad0, key = lambda pt: pt[1]/pt[0])
	
	# Use symmetry to build up the rest of the ordered layer.
	quad1 = [i for i in reversed([(-pt[0],pt[1]) for pt in quad0])]
	semi0 = quad0 + [(0,R)] + quad1
	semi1 = [i for i in reversed([(pt[0],-pt[1]) for pt in semi0])]
	circle = [(R,0)] + semi0 + [(-R,0)] + semi1
	
	# Shift back via origin
	ordered_layer = [(pt[0] + x0, pt[1] + y0) for pt in circle]
	
	return ordered_layer
	
# Return the id's denoting the starts of bins.
def create_bin_partition(layers):
	
	id_prev = 0
	bins = []
	for id in layers:
		num_in_layer = id - id_prev
		numBin = int(math.ceil(num_in_layer*1.0/max_add))
		binSize = int(math.ceil(num_in_layer*1.0/numBin))
		
		bins = bins + [id_prev + j*binSize for j in xrange(1,numBin)]
		bins.append(id_prev + num_in_layer)
		
		id_prev = id
		
	return bins

#- this is somehow necessary, because of some glitch or something.....
# L=4
# H=4

# def compile():
# 	renyiStr = "-DUSE_S:STRING="
# 	gridStr  = "-DUSE_S:STRING=" 
# 	compStr = "cd " + buildDir + "; cmake ../ " + renyiStr + str(Renyi) + " " + gridStr + str(Lattice) + "; make sim;"
# 	print compStr   
# 	subprocess.call( compStr, shell=True)

def main():
	#- Compile
#	compile()

	#- Create the list of directories (just numbers, not actual directories yet)
	if os.path.exists("./"+runDir):
		print "Directory " + runDir + " already exists!"
		print "Exiting script"
		exit(1)
	os.makedirs(runDir)
	print "mkdir " + runDir	
	
	# Loop over the size (L) directories.
	for L in xrange(min_L,max_L,step_L):
		Ldir = "L" + str(L).zfill(3)
		os.makedirs(runDir + "/" + Ldir)
		print "mkdir " + runDir + "/" + Ldir
		
		# Create region data.
		region = RegionData(L,create_layers,create_bin_partition)
		
		#- Loop over the swap ratios (r)
		for r in xrange(0,len(region.bins)):
			# Setup directories
			rDir = "r" + str(r).zfill(3)
			currentDir = runDir+"/" + Ldir + "/" + rDir 
			os.makedirs(currentDir)
			subprocess.call( "cp "+buildDir+"/examples/sim ./"+currentDir+"/"+Ldir+"_"+rDir+".out", shell=True)
			
			#- Create bash_in.txt file & contents
			fname = "%s/bash_in.txt" % currentDir
			f = open(fname, 'w')
			f.write("-shift shift.txt -mult "+str(mult)+" -L "+str(L)+" -H "+str(L)+" -g "+str(r)+"\n")
			f.close()

			#- Create shift.txt file & contents
			sname = "%s/shift.txt" % currentDir
			
			shift = region.generate_shiftmap(L)
# 			printRegion(shift)
			printMatToFile(shift,sname,'w')
			region.set_bin(r,1) # increment by a bin.
			shift = region.generate_shiftmap(L)
# 			printRegion(shift)
			printMatToFile(shift,sname,'a')

main()