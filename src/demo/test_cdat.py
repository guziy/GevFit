__author__="huziy"
__date__ ="$Mar 25, 2011 12:03:32 AM$"


#check whether lmoments package is available
import lmoments

if __name__ == "__main__":
    x = [2, 5, 6, 2, 3, 11, 1]
    x.sort()
    print(x)
    x.remove(2)
    print(x)
    print("Hello World")
