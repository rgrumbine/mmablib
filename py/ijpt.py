"""
#working with ij points, i.e., indices to arrays
#  warning: this are fortran-sense indices, i is most rapidly varying index
#Robert Grumbine
#1 June 2018
"""
import sys

def midpoint(x, y, z):
  """
  find the midpoint between two points in ij space
  """
  z.i = (x.i + y.i)/2.0
  z.j = (x.j + y.j)/2.0
  return z


class ijpt:
  """
  Class ijpt -- like mmablib C++ class ijpt
  """

  def __init__(self,i = 0, j = 0):
    self.i = i
    self.j = j

  def show(self, fout = sys.stderr):
    ''' ijpt.show(fout = xxx) -- print out ijpt to optional file xxx '''
    print (self.i,self.j, file=fout)

  def im(self):
    ''' ijpt.im -- return i-1,j '''
    x = ijpt()
    x.i = self.i - 1
    x.j = self.j
    return x
  def ip(self):
    ''' ijpt.ip -- return i+1,j '''
    x = ijpt()
    x.i = self.i + 1
    x.j = self.j
    return x
  def jm(self):
    ''' ijpt.jm -- return i,j-1 '''
    x = ijpt()
    x.i = self.i
    x.j = self.j - 1
    return x
  def jp(self):
    ''' ijpt.jp -- return i,j+1 '''
    x = ijpt()
    x.i = self.i
    x.j = self.j + 1
    return x
