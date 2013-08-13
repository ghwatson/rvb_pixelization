'''
Created on Aug 13, 2013

@author: grant
'''

import math
from bresenham import raster_disk

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
# Notes: The diagonal values get double-counted.  There is some redundancy in
#        the logic in that raster_disks generates the whole disk, but the
#        function below only needs a quadrant.
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