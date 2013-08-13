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

import os
import subprocess
import shiftlogic

# Prints matrix so that x,y is oriented in the conventional manner.
# Looks like:
# y
# ^
# |
# ---> x
# TODO: Not true! But not a problem at the moment.
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

# A data structure for use in filling out shift files in a general manner.
# The user supplies the logic for providing layers, bins and points.
# CONSTRUCTION:
# points: a 1d list of tuples representing the points, in order of shift
#		  value assignment.
# layer_idx: a list of indexes of points defining the layers.
# 		Example: If we have  [p0 | p1 p2 p3 | p4 p5 | p6], where | defines a
# 			 	 layer, then layer_idx = [1,4,6,7].
# bin_idx: Similar to layer_idx, but instead the indices defining the bins
# 		 used for a single swap ratio.
# shift_state: A list of the shift values of the points.
#
# NOTE: This has only been tested for the disk layering. Test before using
# in other ways.
class RegionData(object):
	def __init__(self,L,max_add,create_layers,create_bin_partition):
		# Get the layer_idx (ex shells of disk, or columns of square).
		self.layer_idx, self.points = create_layers(L)
		self.bin_idx = create_bin_partition(self.layer_idx,max_add)
		
		# Initialize the shift shift_state of the points to 0.
		self.shift_state = [0 for x in xrange(len(self.points))]
	
	# Get the rth bin of shift_state.
	def get_bin(self, r):
		if r is not 0:
			begin = self.bin_idx[r-1]
			end = self.bin_idx[r]
		else:
			begin = end = 0
		return self.shift_state[begin:end]
	
	# Set the rth bin to value.
	def set_bin(self, r, value):
		if r is not 0:
			begin = self.bin_idx[r-1]
			end = self.bin_idx[r]
		else:
			begin = 0
			end = 1
		self.shift_state[begin:end] = [value for x in xrange(end-begin)]
		
	# Returns a 2darray of 0s and 1s that represents the current cartesian 
	# state of the data.
	def generate_shiftmap(self,canvas_L):
		canvas = [[0 for x in xrange(canvas_L)] for x in xrange(canvas_L)]
		for i, value in enumerate(self.shift_state):
			# Break once a 0 is reached.
			if not value:
				break
			# Set the corresponding point of value to 1.
			(x,y) = self.points[i]
			canvas[x][y] = value
		return canvas

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
		region = RegionData(L,max_add,shiftlogic.create_layers,shiftlogic.create_bin_partition)
		
		#- Loop over the swap ratios (r)
		for r in xrange(0,len(region.bin_idx)):
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
#  			printRegion(shift)
			printMatToFile(shift,sname,'w')
			region.set_bin(r,1) # increment by a bin.
			shift = region.generate_shiftmap(L)
#  			printRegion(shift)
			printMatToFile(shift,sname,'a')

main()