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