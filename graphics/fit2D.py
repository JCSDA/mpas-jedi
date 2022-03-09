from copy import deepcopy
import numpy as np
from typing import Iterable, Callable, Union, Tuple, Type, Mapping
#import yaml, sys

Real = Union[float, np.float64, np.single]
Numerical = Union[int, Real]

rVector = Iterable[Real]
coordDict = Mapping[str, Iterable[Numerical]]
termsTuple = Tuple[Iterable[coordDict], Iterable[str]]
termsVectors = Tuple[Iterable[str], rVector]
termsDict = Mapping[str, Iterable[dict]]

## base classes for 2d fitting
class fit2DTerm():
  '''
  single term in a fit2DEquation
  '''
  def __init__(self):
    '''virtual method, initialize constants used in function of x and y'''
    raise NotImplementedError()

  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    '''virtual method, function of x and y'''
    raise NotImplementedError()

  def tname(self) -> Tuple[coordDict, str]:
    '''virtual method, returns tuple like (2-D coordinate tuple, string descriptor of the function)'''
    raise NotImplementedError()


class fit2DEquation():
  '''
  contains multiple terms to be combined by addition or subtraction only in a fit2D object
  '''
  def __init__():
    self._terms = np.asarray([], dtype=fit2DTerm)

  def terms(self) -> Iterable[fit2DTerm]:
    return self._terms

  def termValues(self,
      x: rVector,
      y: rVector,
    ) -> Iterable[rVector]:

    '''produces an array of all self._terms, evaluated at (x, y)'''
    out = []
    for term in self._terms:
      assert isinstance(term, fit2DTerm), 'fit2DEquation.termValues: term must be fit2DTerm'
      out.append(term.evaluate(x, y))
    return np.array(out)

  def termNames(self) -> termsTuple:
    '''produces a list of str's representing the names of all terms'''
    coords = []
    names = []
    for term in self._terms:
      assert isinstance(term, fit2DTerm), 'fit2DEquation.termNames: term must be fit2DTerm'
      coord, n = term.tname()
      coords.append(deepcopy(coord))
      names.append(n)

    return coords, names


### derived classes for 2d fitting
## polynomial base classes
class twoConstant2DTermBase(fit2DTerm):
  '''
  single term in a polynomial-like equation
  Examples:
    identity functions: x^a * y^b
    functions of x^a, y^b: f(x^a)g(y^b)
  '''
  constantNames = 'constants'
  def __init__(self,
      a: Numerical,
      b: Numerical,
    ):
    '''
    Args:
      a: exponent applied to x
      b: exponent applied to y
    '''
    self._a = a
    self._b = b
    self._coord = {self.constantNames: list([int(a), int(b)])}
#    print('test')
#    yaml.dump(
#      self._coord, sys.stdout,
#      #{'terms': [{'exponents': [int(a), int(b)]}]}, sys.stdout,
#      #indent=4,
#      #block_seq_indent=2,
#      allow_unicode=False,
#      #default_flow_style=False,
#    )


class twoConstant2DOrder(fit2DEquation):
  '''
  2D polynomial equation containing all terms of `order`-order only
  Examples for order==4:
    polynomial functions: x^4, x^3y, x^2y^2, xy^3, y^4
    functions of x, y, and algebraicConstants: f(x,4), f(x,3)g(y), f(x,2)g(y,2), f(x)g(y,3), g(y,4)
  '''
  def __init__(self, order: int, termClass: Type[twoConstant2DTermBase]):
    xorders = np.flip(np.arange(0, order+1))
    self._terms = np.empty(len(xorders), dtype=fit2DTerm)
    for cx in xorders:
      cy = order-cx
      self._terms[cx] = termClass(float(cx), float(cy))


class twoConstant2DEquationBase(fit2DEquation):
  '''
  2D polynomial equation containing all terms of maxorder and smaller
  '''
  termClass = twoConstant2DTermBase
  def __init__(self, maxorder: int=2):
     orders = np.arange(0, maxorder+1)
     self._terms = np.asarray([], dtype=fit2DTerm)
     for o in orders:
       orderTerms = twoConstant2DOrder(o, self.termClass)
       self._terms = np.append(self._terms, orderTerms.terms())


## specific polynomial families
# classic polynomial
class poly2DTerm(twoConstant2DTermBase):
  '''evaluates x^a * y^b'''
  constantNames = 'exponents'
  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    return np.power(x, self._a) * np.power(y, self._b)

  def tname(self) -> Tuple[coordDict, str]:
    out = ''
    cx = self._a
    cy = self._b
    if cx==0 and cy == 0:
      out = 'c'
    else:
      if cx > 1:
        out += 'x^'+str(cx)
      elif cx == 1:
        out += 'x'
      if cy > 1:
        out += 'y^'+str(cy)
      elif cy == 1:
        out += 'y'

    return self._coord, out


class poly2DEquation(twoConstant2DEquationBase):
  '''twoConstant2DEquationBase derived class for a poly2DTerm'''
  termClass = poly2DTerm


# products of (1-polynomial)
class polyOneSub2DTerm(twoConstant2DTermBase):
  '''evaluates (1 - x^a) * (1 - y^b)'''
  constantNames = 'exponents'
  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    return (1. - np.power(x, self._a)) * (1. - np.power(y, self._b))

  def tname(self) -> Tuple[coordDict, str]:
    out = ''
    cx = self._a
    cy = self._b
    if cx==0 and cy == 0:
      out = 'c'
    else:
      if cx > 1:
        out += '(1-x^'+str(cx)+')'
      elif cx == 1:
        out += '(1-x)'
      if cy > 1:
        out += '(1-y^'+str(cy)+')'
      elif cy == 1:
        out += '(1-y)'

    return self._coord, out


class polyOneSub2DEquation(twoConstant2DEquationBase):
  '''twoConstant2DEquationBase derived class for a polyOneSub2DTerm'''
  termClass = polyOneSub2DTerm


## specific trigonometric families
# not very good
class cosineHalfPiProduct2DTerm(twoConstant2DTermBase):
  '''evaluates cos(a * pi/2 * x) * cos(b * pi/2 * y)'''
  constantNames = 'coefficients of dependent variables'
  def __init__(self,
      a: Numerical,
      b: Numerical,
    ):
    super().__init__(a, b)
    self.__px = float(self._a) * np.pi * 0.5
    self.__py = float(self._b) * np.pi * 0.5

  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    return np.cos(self.__px * x) * np.cos(self.__py * y)

  def tname(self) -> Tuple[coordDict, str]:
    out = ''
    cx = self._a
    cy = self._b
    if cx==0 and cy == 0:
      out = 'c'
    else:
      if cx > 1:
        out += 'cos('+str(cx)+'*pi/2*x)'
      elif cx == 1:
        out += 'cos(pi/2*x)'
      if cy > 1:
        out += 'cos('+str(cy)+'*pi/2*y)'
      elif cy == 1:
        out += 'cos(pi/2*y)'

    return self._coord, out


class cosineHalfPiProduct2DEquation(twoConstant2DEquationBase):
  '''twoConstant2DEquationBase derived class for a cosineHalfPiProduct2DTerm'''
  termClass = cosineHalfPiProduct2DTerm


# also not very good
class sineHalfPiPoly2DTerm(twoConstant2DTermBase):
  '''evaluates sin(pi/2 * x^a) * sin(pi/2 * y^b)'''
  constantNames = 'exponents'
  halfpi = np.pi * 0.5
  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    return np.sin(self.halfpi * np.power(x, self._a)) * np.sin(self.halfpi * np.power(y, self._b))

  def tname(self) -> Tuple[coordDict, str]:
    out = ''
    cx = self._a
    cy = self._b
    if cx==0 and cy == 0:
      out = 'c'
    else:
      if cx > 1:
        out += 'sin(pi/2*x^'+str(cx)+')'
      elif cx == 1:
        out += 'sin(pi/2*x)'
      if cy > 1:
        out += 'sin(pi/2*y^'+str(cy)+')'
      elif cy == 1:
        out += 'sin(pi/2*y)'

    return self._coord, out


class sineHalfPiPoly2DEquation(twoConstant2DEquationBase):
  '''twoConstant2DEquationBase derived class for a sineHalfPiPoly2DTerm'''
  termClass = sineHalfPiPoly2DTerm


# did not converge in least squares for CFy,CFx problem
class expBinomial2DTerm(twoConstant2DTermBase):
  '''evaluates exp(x^a) * exp(y^b) == exp(x^a + y^b)'''
  constantNames = 'exponents'
  def evaluate(self,
      x: rVector,
      y: rVector,
    ) -> rVector:
    return np.exp(np.power(x, self._a) + np.power(y, self._b))

  def tname(self) -> Tuple[coordDict, str]:
    out = ''
    cx = self._a
    cy = self._b
    if cx==0 and cy == 0:
      out = 'c'
    else:
      if cx > 1:
        #out += 'exp('+str(cx)+'*x)'
        out += 'exp(x^'+str(cx)+')'
      elif cx == 1:
        out += 'exp(x)'
      if cy > 1:
        #out += 'exp('+str(cy)+'*y)'
        out += 'exp(y^'+str(cy)+')'
      elif cy == 1:
        out += 'exp(y)'

    return self._coord, out


class expBinomial2DEquation(twoConstant2DEquationBase):
  '''twoConstant2DEquationBase derived class for a expBinomial2DTerm'''
  termClass = expBinomial2DTerm


## 2D fitting class
class fit2D():
  '''
  utility for fitting a general 2d equation to 3-dimensional data
  defaults to a degree-2 2d polynomial
  '''
  def __init__(self,
      xd: rVector,
      yd: rVector,
      zd: rVector,
      equation: fit2DEquation = poly2DEquation(2),
      xytransform: Callable[[rVector], rVector]=None,
    ):
    '''
    Args:
      xd: x predictor data
      yd: y predictor data
      zd: data to fit
      equation: equation with which to fit z-data
      xytransform: Callable one-to-one transformation to be applied to x and y before
                   fitting to and predicting z data
    '''
    if xytransform is None:
      self.__transform = self.identityTransform
    elif callable(xytransform):
      self.__transform = xytransform
    else:
      print('WARNING in fit2D.__init__: xytransform is not callable, using self.identityTransform instead')
      self.__transform = self.identityTransform

    xt = self.__transform(xd)
    yt = self.__transform(yd)
    self.__equation = equation
    terms = self.__equation.termValues(xt, yt)
    self.__nterms = len(terms)

    c, r, rank, s = np.linalg.lstsq(terms.T, zd, rcond=None)
    self.__coeffs = np.asarray(c, dtype=np.float64)

  @staticmethod
  def identityTransform(x: rVector) -> rVector:
    return x

  def predict(self,
      x: rVector,
      y: rVector,
    ) -> rVector:

    '''evaluates the 2D polynomial at (x, y)'''
    #out = np.zeros(x.shape)
    xt = self.__transform(np.asarray(x).flatten())
    yt = self.__transform(np.asarray(y).flatten())
    terms = self.__equation.termValues(xt, yt)

    # create intermediate result with terms multiplied by coefficients
    fullTerms = np.empty((self.__nterms, len(xt)))
    for iterm, c in enumerate(self.__coeffs):
      fullTerms[iterm, :] = c * terms[iterm]

    # sum over first axis, return to starting shape
    out = fullTerms.sum(axis=0).reshape(np.asarray(x).shape)

    return out

  def coeffs(self):
    '''returns the polynomial coefficients for this instance'''
    return deepcopy(self.__coeffs)

  def equation(self, precision: int=3):
    '''produces a dictionary of formatted fit2DEquation terms'''
    out = {}
    #fmt = ' {:0.'+str(precision)+'f}'
    fmt = ' {:0.'+str(precision)+'f}'
    coords, names = self.__equation.termNames()
    for ic, (coord, name) in enumerate(zip(coords, names)):
        c = self.__coeffs[ic]
        out[coord] = {
          'cstr':'',
          'c': c,
        }
        out[coord]['name'] = str(name)

        if ic > 0 and self.__coeffs[ic] > 0: out[coord]['cstr']+='+'
        out[coord]['cstr'] = out[coord]['cstr'] + fmt.format(c)
        out[coord]['cstr'] = out[coord]['cstr'].replace(' ', '')

    return out

  def terms(self, precision: int=3, sortby: str=None, returnType='list') -> Union[termsVectors, termsDict]:
    fmt = ' {:0.'+str(precision)+'e}'
    coords, names = self.__equation.termNames()
    cVals = np.asarray(self.coeffs())
    cFormatted = np.asarray([fmt.format(c).replace(' ', '') for c in cVals])

    coeffSort = np.arange(len(cVals))
    if sortby is not None:
      if sortby == 'coefficient':
        coeffSort = np.flip(np.argsort(np.abs(cVals)))
      elif sortby == 'names' and isinstance(names[0], str):
        coeffSort = np.flip(np.argsort(np.asarray(names)))

    if returnType == 'list':
      return np.asarray(names)[coeffSort], cFormatted[coeffSort]

    elif returnType == 'dict':
      eqDict = {'terms': []}
      for ic in coeffSort:
        coord = coords[ic]
        if isinstance(coord, dict) and len(coord) > 0:
          d = deepcopy(coord)
          d['value'] = float(cVals[ic])
        else:
          d = {'name': str(names[ic]), 'value': float(cVals[ic])}
        eqDict['terms'].append(d)
      return eqDict

    else:
      return None, None
