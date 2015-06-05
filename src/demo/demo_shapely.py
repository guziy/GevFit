# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="huziy"
__date__ ="$Aug 10, 2011 8:59:27 AM$"

from shapely.geometry import Polygon, Point


if __name__ == "__main__":
    p = Point(5,5)
    polygon = Polygon()

    print p.to_wkt()
    print "Hello World"
