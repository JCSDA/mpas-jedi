#!/usr/bin/env python3

import itertools as itt
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np

#============
# plot styles
#============

## cvdCompliant: whether color palettes should be color vision deficiency (CVD) compliant
cvdCompliant = True

## lowFidelity: whether to use coarser color map bins for 2D images
lowFidelity = False

## styleOrder
# options:
# + 'both' [both color and line iterate simultaneously,
#           as long as their respective iterator lengths
#           are mutually non-divisible]
#    example:
#      colorIterator = ['b', 'g', 'm'], lineIterator = ['-', '--', '-.', ':']
#    produces:
#      orderedColors = ['b', 'g', 'm', 'b', ..., 'm', 'b', 'g', 'm']
#      orderedLines = ['-', '--', '-.', ':', ..., '-', '--', '-.', ':']
#
# + 'both-force' [both color and line iterate simultaneously,
#                 regardless of their lengths]
#
# + 'colors then lines' [fast varying colors]
#   example:
#     colorIterator = ['b', 'g', 'm'], lineIterator = ['-', '--', '-.', ':']
#   produces:
#     orderedColors = ['b', 'g', 'm', ..., 'b', 'g', 'm']
#     orderedLines = ['-', '-', '-', ..., ':', ':', ':']
#
# + 'lines then colors' [fast varying lines]
#   example:
#     colorIterator = ['b', 'g', 'm'], lineIterator = ['-', '--', '-.', ':']
#   produces:
#     orderedColors = ['b', 'b', 'b', 'b', ..., 'm', 'm', 'm', 'm']
#     orderedLines = ['-', '--', '-.', ':', ..., '-', '--', '-.', ':']
#
# + otherwise, default is equivalent to 'colorIterator[0] and all lines'
styleOrder = 'both'

spaghettiColor = ['0.45']*3
spaghettiLine = '--'

# nSpaghetti:
# forces the leading nSpaghetti styles to take on the spaghetti characteristics
# default(0) is no spaghetti
nSpaghetti = 0
assert nSpaghetti >= 0, 'nSpaghetti is positive semi-definite'

## line styles
lineIterator = ['-', '--', '-.', ':']

if cvdCompliant:
  ## 5-color, custom-generated using multiple palettes from https://colorbrewer2.org
  # checked for CVD-compliance at https://www.color-blindness.com/coblis-color-blindness-simulator/
  # black added after the fact, uncomment for CONTROL experiment
  #colorIterator = [
  #  #'#000',
  #  '#66c2a5',
  #  '#1f78b4',
  #  '#b2df8a',
  #  '#1b9e77',
  #  '#d95f02',
  #]

  ## 3-color, compliant for color vision deficiency (CVD) generated at https://colorbrewer2.org
  # black added after the fact, uncomment for CONTROL experiment
  colorIterator = [
  #  #'#000',
    '#d95f02', # dark orange
    '#7570b3', # purple
    '#1b9e77', # green
  ]

  ## 4-color, compliant for color vision deficiency (CVD) generated at https://colorbrewer2.org
  # black added after the fact, uncomment for CONTROL experiment
  #colorIterator = [
  ##  '#000',
  #  '#a6cee3', # light blue
  #  '#1f78b4', # blue
  #  '#b2df8a', # light green
  #  '#33a02c', # green
  #]
else:
  # 7-color, not CVD-compliant
  colorIterator = [
    [0., 0., 0.],
    [0.0000, 0.4470, 0.7410], # blue
    [0.8500, 0.3250, 0.0980], # red
    [0.9290, 0.6940, 0.1250], # gold
    [0.4940, 0.1840, 0.5560], # purple
    [0.4660, 0.6740, 0.1880], # green
    [0.3010, 0.7450, 0.9330], # light blue
    [0.6350, 0.0780, 0.1840], # crimson
  ]

## grey
#colorIterator = [
#  [0.2, 0.2, 0.2],
#  [0.4, 0.4, 0.4],
#  [0.6, 0.6, 0.6],
#  [0.8, 0.8, 0.8],
#]

nColor = len(colorIterator)
nLine = len(lineIterator)

if (
  (styleOrder == 'both' and
  np.mod(nColor, nLine) > 0 and
  np.mod(nLine, nColor) > 0)
  or styleOrder == 'both-force'
  ):
  # works best when line and color iterator lengths are mutually non-divisible
  orderedColors = colorIterator
  orderedLines = lineIterator
else:
  orderedLines = []
  orderedColors = []
  if styleOrder == 'colors then lines':
    styleGenerator = itt.product(
      lineIterator,
      colorIterator)
    for (l, c) in styleGenerator:
      orderedColors.append(c)
      orderedLines.append(l)

  elif styleOrder == 'lines then colors':
    styleGenerator = itt.product(
      colorIterator,
      lineIterator)
    for (c, l) in styleGenerator:
      orderedColors.append(c)
      orderedLines.append(l)

  else:
    orderedLines = lineIterator
    orderedColors = [colorIterator[0]]

def nth(iterable, n, default=None):
  "Returns the nth item or a default value"
  return next(itt.islice(iterable, n, None), default)

def color(index: int):
  assert index >=0, 'index is positive semi-definite'
  if index < nSpaghetti:
    return spaghettiColor
  else:
    colorCycle = itt.cycle(orderedColors)
    return nth(colorCycle, index-nSpaghetti)

def line(index: int):
  assert index >=0, 'index is positive semi-definite'
  if index < nSpaghetti:
    return spaghettiLine
  else:
    lineCycle = itt.cycle(orderedLines)
    return nth(lineCycle, index-nSpaghetti)


## color maps for 2D figures
# color maps compliant for color vision deficiency (CVD) generated at https://colorbrewer2.org
# sequential
# must be even number to avoid errors with BoundaryNorm
OrRd10 = ['#fff','#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000']
Reds10 = ['#fff','#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d']
YlOrBr10 = ['#fff','#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506']
PuRd10 = ['#fff','#f7f4f9','#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f']
sequential = PuRd10
sequentialCMap = mplcolors.ListedColormap(sequential)

# diverging
# must be even number to place 0.0 properly at origin of symmetric diverging colors
RdBu10 = ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
RdYlBu10 = ['#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695']
#RdYlBu11 = ['#a50026','#d73027','#f46d43','#fdae61','#fee090','#fff','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695']

PRGn10 = ['#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b']
diverging = RdYlBu10
diverging.reverse() # for warm-to-cold maps to be changed to cold-to-warm (e.g., swap red and blue)
divergingCMap = mplcolors.ListedColormap(diverging)

# color map and # of levels
if lowFidelity:
  cmaps = {
    'sequential': {
      'map': sequentialCMap,
      'n': sequentialCMap.N,
     },
    'diverging': {
      'map': divergingCMap,
      'n': divergingCMap.N,
     },
  }
else:
  cmaps = {
    'sequential': {
      'map': plt.get_cmap('YlOrBr'),
      #'map': plt.get_cmap('purd'),
      #'map': plt.get_cmap('YlGnBu'),
      'n': 18,
     },
    'diverging': {
      #'map': plt.get_cmap('seismic'),
      'map': plt.get_cmap('bwr'),
      'n': 26,
     },
  }


## subplot title labels
subplotIterator = [
'a', 'b', 'c', 'd', 'e',
'f', 'g', 'h', 'i', 'j',
'k', 'l', 'm', 'n', 'o',
'p', 'q', 'r', 's', 't',
'y', 'v', 'w', 'x', 'y',
'z',
'A', 'B', 'C', 'D', 'E',
'F', 'G', 'H', 'I', 'J',
'K', 'L', 'M', 'N', 'O',
'P', 'Q', 'R', 'S', 'T',
'Y', 'V', 'W', 'X', 'Y',
'Z',
'aa', 'ab', 'ac', 'ad', 'ae',
'af', 'ag', 'ah', 'ai', 'aj',
'ak', 'al', 'am', 'an', 'ao',
'ap', 'aq', 'ar', 'as', 'at',
'ay', 'av', 'aw', 'ax', 'ay',
'az',
'ba', 'bb', 'bc', 'bd', 'be',
'bf', 'bg', 'bh', 'bi', 'bj',
'bk', 'bl', 'bm', 'bn', 'bo',
'bp', 'bq', 'br', 'bs', 'bt',
'by', 'bv', 'bw', 'bx', 'by',
'bz',
'ca', 'cb', 'cc', 'cd', 'ce',
'cf', 'cg', 'ch', 'ci', 'cj',
'ck', 'cl', 'cm', 'cn', 'co',
'cp', 'cq', 'cr', 'cs', 'ct',
'cy', 'cv', 'cw', 'cx', 'cy',
'cz',
'da', 'db', 'dc', 'dd', 'de',
'df', 'dg', 'dh', 'di', 'dj',
'dk', 'dl', 'dm', 'dn', 'do',
'dp', 'dq', 'dr', 'ds', 'dt',
'dy', 'dv', 'dw', 'dx', 'dy',
'dz',
]

def subplotLabel(index = 0):
  if index < 0:
    return ''
  else:
    return '('+subplotIterator[np.mod(index, len(subplotIterator))]+')'


## other line and color styles
plotSpecs = ['k-*', 'b-*', 'g-*', 'r-*', 'c-*', 'm-*',
             'k--+','b--+','g--+','r--+','c--+','m--+']
plotMarkers = ['*','*','*','*','*','*',
               '+','+','+','+','+','+']
