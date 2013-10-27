'''
Created on Aug 13, 2013

Provides the logic for filling out disks for use in grantscript.

This is essentially the same as shiftlogic.py, but instead does for
off-site disks with odd numbers of pixels in the diameter (i.e. the radius
and the origin are half-integer). So, there are a variety of tweaks.

Note: As of now, I just hacked the original shiftlogic.py to do this.
There might be a better way to implement this stuff. Ie, it'd be ideal
to just have one shiftlogic.py that works in general for on-site or off-site
circles.

@author: grant
'''

import math
from bresenham import raster_disk_offsite_integer_rad

#removes duplicates in a list whilst preserving order.
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
	
	x0 = x0 - 0.5 # Shift to the upper right offsite point.
	y0 = y0 - 0.5
	max_r = max_r + 0.5 # We can fit an extra 0.5.
	
	# Collect all the disks full of points.
	rasterdisks = {}
	for d in xrange(1,int(2*max_r),2):
		rasterdisks[d] = raster_disk_offsite_integer_rad(x0,y0,d*1./2)

	# Build the layers.
	for d in xrange(1,int(2*max_r),2):
		#TODO: don't hardcode the step in the diameters. ugly hack.
		if d > 2:
			layer_of_points = [x for x in rasterdisks[d] if x not in rasterdisks[d-2]]
		else:
			layer_of_points = rasterdisks[d]
		layer_of_points = remove_duplicates(layer_of_points)
		layer_of_points = order_disk_layer(layer_of_points,d*1./2,(x0,y0))
		layer_of_points = remove_duplicates(layer_of_points)
		
		points = points + layer_of_points
		layer_idx.append(len(points))
		
	return layer_idx, points

# Order the points in a layer starting from (R,0) and going around CCW.
# Notes: The diagonal values get double-counted.  There is some redundancy in
#        the logic in that raster_disks generates the whole disk, but the
#        function below only needs a quadrant.
# TODO: R is not needed!
def order_disk_layer(layer,R,origin):
	(x0,y0) = origin
	# Move points to origin
	layer = [(pt[0] - x0,pt[1] - y0) for pt in layer]
	
	# Remove points not in first quadrant.
	quad0 = [ pt for pt in layer if (pt[0] > 0 and pt[1] > 0) ]
		
	# Perform sorting.
	quad0 = sorted(quad0, key = lambda pt: pt[1]/pt[0])
	
	# Use symmetry to build up the rest of the ordered layer.
	quad1 = [i for i in reversed([(-pt[0],pt[1]) for pt in quad0])]
	semi0 = quad0 + quad1
	semi1 = [i for i in reversed([(pt[0],-pt[1]) for pt in semi0])]
	circle = semi0 + semi1
	
	# Shift back via origin
	ordered_layer = [(pt[0] + x0, pt[1] + y0) for pt in circle]
	#Again, ugly anti-pythonic int conversion hacking.
	#TODO: fix this
	ordered_layer = [(int(x),int(y)) for (x,y) in ordered_layer]
	
	return ordered_layer
	
# Return the id's denoting the starts of bins.
def create_bin_partition(layers, max_add):
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