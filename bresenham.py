'''
Created on Aug 13, 2013

A module providing Bresenham's circle algorithm, altered to provide a disk.
Ref: http://en.wikipedia.org/wiki/Midpoint_circle_algorithm

@author: grant
'''

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
#     http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
# The ordering of the points in the set is like so:
#          13 12 11   
#       10  9  8  7  6
#        5  4  3  2  1
#       18 17 16 15 14
#          21 20 19
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

# This algorithm is capable of dealing with half-integer origin and integer
# radii, unlike raster_disk which can only deal with integer origin and integer
# radii.
# TODO: 
# -Get rid of the int conversion hack by redesigning stuff as needed.
def raster_disk_offsite_integer_rad(x0,y0,radius):
    x = 0.5 # Start with the rightmost of the 2 upper pixels.
    y = radius
    
    coords = []
    
    f = 1 - radius
    ddF_x = 1
    ddF_y = -2 * radius
    
    
    while (y >= x):
        #TODO: This conversion is pretty anti-pythonic.
        x1 = int(x0+x)
        x2 = int(x0-x)
        x3 = int(x0 + y)
        x4 = int(x0 - y)
        y1 = int(y0 + y)
        y2 = int(y0 - y)
        y3 = int(y0 + x)
        y4 = int(y0 - x)
        
        line1 = get_horizontal_line((x1, y1), (x2, y1))
        line2 = get_horizontal_line((x1, y2), (x2, y2))
        line3 = get_horizontal_line((x3, y3), (x4, y3))
        line4 = get_horizontal_line((x3, y4), (x4, y4))
         
        coords = coords + line1 + line2 + line3 + line4
        
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
        
    return coords