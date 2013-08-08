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
#parameter["sq"]    = "-r 30h --mpp 1g"
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

# Generate the coordinates for a circle at pixel (x0,y0) with integer radius.
# Script to get pixelized circle derived from wikipedia at:
#     http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
def raster_circle(x0, y0, radius):
    # Generate the coordinates.
    f = 1 - radius
    ddF_x = 1
    ddF_y = -2 * radius
    x = 0
    y = radius
    
    coords = []
    
    coords.append((x0, y0 + radius))
    coords.append((x0, y0 - radius))
    coords.append((x0 + radius, y0))
    coords.append((x0 - radius, y0))
    
    while x < y :
#         ddF_x == 2 * x + 1
#         ddF_y == -2 * y
#         f == x*x + y*y - radius*radius + 2*x - y + 1
        if f >= 0:
            y -= 1
            ddF_y += 2
            f += ddF_y
        x += 1
        ddF_x += 2
        f += ddF_x    

        coords.append((x0 + x, y0 + y))
        coords.append((x0 - x, y0 + y))
        coords.append((x0 + x, y0 - y))
        coords.append((x0 - x, y0 - y))
        coords.append((x0 + y, y0 + x))
        coords.append((x0 - y, y0 + x))
        coords.append((x0 + y, y0 - x))
        coords.append((x0 - y, y0 - x))
    
    return coords

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
# of the Bresenham circle algorithm.
# The ordering of the points in the set is like so:
#              9
#           8  7  6
#        5  4  3  2  1
#          10 11  12
#             13

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
#         ddF_x == 2 * x + 1
#         ddF_y == -2 * y
#         f == x*x + y*y - radius*radius + 2*x - y + 1
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

#- generate the correct number of folders for an L/2 by L/2 square region 
#- with the given maximum number of sites added per ratio-step
def Bins(L_):
	R_max = L_/4
	numBin =  int(math.ceil( (L_/2)*1.0/max_add ))
	allbins = [l for l in range(0, numBin*L_/2)]
	return allbins

# Add up the bins of all shells making up the largest circle with the 
# L/2 x L/2 as the bounding box.
def TotalBins(L_):
	r_bound = int(math.floor(L_/4))
	numBins = 0
	for r in range(0,r_bound):
		numBins += BinsOfShell(r, L_)
	return numBins

def BinsOfShell(r, L_):
	# Figure out the number of bins needed just for this shell.
	# 0. Define origin as middle of the L/2 box.
	# If even number of nodes, then we opt for slightly smaller box, or
	# slighly larger? For now, lets say the smaller box.
	# 1. Draw out the circle, given the origin.
	# Drawing all these circles will probably be needed later...?
	# These will give the pixel coordinates making up each circle, and thus
	# the path to add along.
	# Perhaps calculate this at the start and feed into TotalBins.
	# 2. Perform some of Ann's stuff.
	x0 = y0 = int(math.floor(L_/4))
	coords = rasterCircle(x0,y0,r)
	
	numBins = int(math.ceil( len(coords)*1.0/max_add))
	return numBins

# Generate the coordinates of the disk, in concentric circles.
# (what to do about non-integer radii?)
def generateDisk(R):
	circles = {}
	for r in range(0,R):
		x0 = y0 = int(math.floor(L_/4))
		circles[r] = rasterCircle(x0,y0,r)

# Gives you the vector describing the bins to fill a row.
def addVector(L_,max_add_):
	numBin =  int(math.ceil( (L_/2)*1.0/max_add ))
	binSize = int(math.ceil((L_/2)*1.0/numBin))
	vAdd = [binSize]*(numBin-1) 
	if ((L_/2)%binSize)>0:
		vAdd.append((L_/2)%(binSize*(numBin-1)))
	else:
		vAdd.append(binSize)
	return vAdd

def addVectorCircle(r,max_add_):
	pass

#- Add sites to region A of the swap or preswap
def regionAadd(L_,r_,addVec):
	tally = 0
	for i in range(0,r_):
		a = i%len(addVec)
		tally += addVec[a]
	Swap = ['1']*tally
	for i in range(0,L_/2*L_/2-tally):
		Swap.append('0')
	fullSwap = []
	tempSwap = []
	for i in range(0,L_/2):
		for j in range(0,L_/2):
			tempSwap.append(Swap[i*L_/2+j])
		for j in range(0,L_/2):
			tempSwap.append('0')
		fullSwap.append(tempSwap)
		tempSwap=[]
	for i in range(0,L_/2):
		for i in range(0,L_):
			tempSwap.append('0')
		fullSwap.append(tempSwap)
		tempSwap=[]
	#printMatrix(fullSwap)
	return fullSwap

def circleAadd(L_,r_,addVec):
	# get a tally of how many 1s there are already.
	tally = 0
	for i in range(0,r_):
		a = i%len(addVec)
		tally += addVec[a]
	
	# Get the full 1d profile of the region A	
	Swap = ['1']*tally
	for i in range(0,L_/2*L_/2-tally):
		Swap.append('0')
	
	fullSwap = []
	tempSwap = []
	for i in range(0,L_/2):
		for j in range(0,L_/2):
			tempSwap.append(Swap[i*L_/2+j])
		for j in range(0,L_/2):
			tempSwap.append('0')
		fullSwap.append(tempSwap)
		tempSwap=[]
	for i in range(0,L_/2):
		for i in range(0,L_):
			tempSwap.append('0')
		fullSwap.append(tempSwap)
		tempSwap=[]
	#printMatrix(fullSwap)
	return fullSwap
	
def printMatrix(mat):
	for bla in mat:
	    	print ' '.join(map(str, bla))
	print('\n')

def printMatToFile(mat,filename,option):
	m = open(filename, option)
	for i in mat:
		m.write(" ".join(i) + "\n") 
	m.write('\n')
	m.close()
	
# TODO: add some more classes. ex: shift-point, which contains a pt and a
#         shift value.
class LayeredBinsData(object):
	def __init__(self,L,create_layers,create_bin_partition):
		# Get the layers (ex shells of disk, or columns of square).
		self.layers = create_layers(L)
		
		for key in self.layers.keys():            
			# Zip in the shift values, all initialized to 0.
			zeros = [[0] for i in xrange(len(self.layers[key]))]
			self.layers[key] = zip(self.layers[key],zeros)
			# Perform binning.
			self.layers[key] = create_bin_partition(self.layers[key])

#removes duplicates in a list whilst preserving order
def remove_duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

# Create circles with integer radii.
# TODO: edit the function rasterCircle to allocate points in the desired
# 		manner...or perform postprocessing here.
def create_layers(L):
	layers = {}
	x0 = y0 = int(math.floor(L/4))
	max_r = int(math.ceil(L/4))
 	layers[0] = [(x0,y0)]

	rasterdisks = {}
	rasterdisks[0] = layers[0]
	for r in xrange(1,max_r):
		rasterdisks[r] = raster_disk(x0,y0,r)

	for r in xrange(1,max_r):
		layers[r] = [x for x in rasterdisks[r] if x not in rasterdisks[r-1]]
		layers[r] = remove_duplicates(layers[r])
	return layers

# Order the points in a layer starting from (R,0) and going around CCW.
# TODO: we only needed the octant here for ordering...consider just working
# with octants as opposed to the full thing.
def order_disk_layer(layer,R):
	# Remove points not in first octant.
	y_bound = R/math.sqrt(2)
	octant0 = [ pt for pt in layer if (pt[0] > 0 and 0 <= y <= y_bound) ]
	
	# Perform sorting by use of a rotated basis.		
	octant0 = sorted(octant0, key = lambda pt:(1/math*sqrt(2) * (pt[0]-pt[1])) )
	
	# Use symmetry to build up the rest of the ordered layer.
	octant1 = [i for i in reversed([(pt[1],pt[0]) for pt in octant0])]
	quadrant0 = octant0 + octant1
	quadrant1 = [i for i in reversed([-pt[0],pt[1]] for pt in quadrant0)]
	semi0 = quadrant0 + quadrant1
	semi1 = [i for i in reversed([pt[0],-pt[1]] for pt in semi0)]
	
	return (semi0 + semi1)
		
		
def create_bin_partition(layer):
	numBin = int(math.ceil(len(layer)*1.0/max_add))
	binSize = int(math.ceil(len(layer)*1.0/numBin))
	
	binned_layer = []
	for i in [j*binSize for j in xrange(numBin-1)]:
		binned_layer.append(layer[i:(i+binSize)])
	binned_layer.append(layer[binSize*(numBin-1):])
	
	return binned_layer

#- this is somehow necessary, because of some glitch or something.....
# L=4
# H=4

def compile():
	renyiStr = "-DUSE_S:STRING="
        gridStr  = "-DUSE_S:STRING=" 
        compStr = "cd " + buildDir + "; cmake ../ " + renyiStr + str(Renyi) + " " + gridStr + str(Lattice) + "; make sim;"
        print compStr   
        subprocess.call( compStr, shell=True)

def testmain():
    L = 50.0
    
#     mylayer = [1,2,5,8,9,10,12,15,30,28,30,38,51,52,90,28,501,202]
#     out = create_bin_partition(mylayer)
    
    # Use the OG code here to compare.
#     max_add = 8
#     avec = addVector(22,max_add)
#     print avec
    
    data = LayeredBinsData(L,create_layers,create_bin_partition)
    
    # Draw each shell.
    for layer in data.layers.values():
        canvas = [ [0 for i in xrange(int(L))] for j in xrange(int(L)) ]
        # Plot each point.
        for bin in layer:
            for shiftpt in bin:
                pt = shiftpt[0]
                canvas[pt[0]][pt[1]] = 1
        printMatrix(canvas)
        pass
    
            
    pass
    print 'ay yo'

def main2():
	#- Compile
#	compile()

	#- Create the list of directories (just numbers, not actual directories yet)
	dirs  = [[l, Bins(l)] for l in range(min_D, max_D, step_D)]
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
		
		
		# Create data structure.
		# structure
		# 	layers
		#		bins
		current_shift = create_data_structure(L)
		
		#- Loop over the swap ratios (r)
		for r in xrange(0,TotalBins(L)):
			# Setup directories
			rDir = "r" + str(r).zfill(3)
			currentDir = runDir+"/" + Ldir + "/" + rDir 
			os.makedirs(currentDir)
			subprocess.call( "cp "+buildDir+"/examples/sim ./"+currentDir+"/"+Ldir+"_"+rDir+".out", shell=True)
			
			#- Create bash_in.txt file & contents
			fname = "%s/bash_in.txt" % currentDir
			f = open(fname, 'w')
			f.write("-shift shift.txt -mult "+str(mult)+" -L "+L+" -H "+L+" -g "+str(r)+"\n")
			f.close()

			#- Create shift.txt file & contents
			sname = "%s/shift.txt" % currentDir
			
			# Create the add vector for the radius that we're interested in.
			aVec = addVector(Lsize,max_add)
			
			shift = circleAadd(Lsize,r,aVec)
 			printMatToFile(shift,sname,'w')
			shift = circleAadd(Lsize,r+1,aVec)
 			printMatToFile(shift,sname,'a')


# def main():
# 
# 	#- Compile
# #	compile()
# 
# 	#- Create the list of directories (just numbers, not actual directories yet)
# 	dirs  = [[l, Bins(l)] for l in range(min_L, max_L, step_L)]
# 	if os.path.exists("./"+runDir):
# 		print "Directory " + runDir + " already exists!"
# 		print "Exiting script"
# 		exit(1)
# 	os.makedirs(runDir)
# 	print "mkdir " + runDir	
# 		
# 	#- Loop over the size (L) directories
# 	for i in range(0,len(dirs)):
# 		Lsize = dirs[i][0]
# 		L = str(Lsize);
# 		Ldir = "L"+L.zfill(3)
# 		os.makedirs(runDir+"/"+Ldir)	
# 		print "mkdir " + runDir+"/"+Ldir
# 
# 		aVec = addVector(Lsize,max_add)
# 
# 		#- Loop over the swap ratios (r)
# 		for r in dirs[i][1]:
# 			rDir = "r"+str(r).zfill(3)
# 			currentDir = runDir+"/"+Ldir+"/"+rDir 
# 			os.makedirs(currentDir)
# 			#print "mkdir " + runDir+"/"+Ldir+"/"+rDir
# 			subprocess.call( "cp "+buildDir+"/examples/sim ./"+currentDir+"/"+Ldir+"_"+rDir+".out", shell=True)
# 			
# 			#- Create bash_in.txt file & contents
# 			fname = "%s/bash_in.txt" % currentDir
# 			f = open(fname, 'w')
# 			f.write("-shift shift.txt -mult "+str(mult)+" -L "+L+" -H "+L+" -g "+str(r)+"\n")
# 			f.close()
# 
# 			#- Create shift.txt file & contents
# 			sname = "%s/shift.txt" % currentDir
# 			
# 			shift = regionAadd(Lsize,r,aVec)
#  			printMatToFile(shift,sname,'w')
# 			shift = regionAadd(Lsize,r+1,aVec)
#  			printMatToFile(shift,sname,'a')
#  			
# 			shift = circleAadd(Lsize,r,aVec)
#  			printMatToFile(shift,sname,'w')
# 			shift = circleAadd(Lsize,r+1,aVec)
#  			printMatToFile(shift,sname,'a')
# 			




testmain()
