#!/usr/bin/env python3

from collections.abc import Iterable
from copy import deepcopy
import os
import itertools
import logging
import multiprocessing as mp
import numpy as np
from scipy.spatial import cKDTree

_logger = logging.getLogger(__name__)

class InterpolateCartesian:
  eps = 1.e-10
  def __init__(self, XIn, YIn, ZIn, nInterpPoints = 5,
               weightMethod = 'barycentricScalar',
               distanceMethod = 'cartesian',
               calculateDiagnostics = False):
    assert len(XIn) == len(YIn), 'InterpolateCartesian: input coordinates must be identical size'
    assert len(XIn) == len(ZIn), 'InterpolateCartesian: input coordinates must be identical size'
    self.calculateDiagnostics = calculateDiagnostics

    # init weight method
    weightMethods = {
      'unstinterpScalar': {
        'method': self.unstinterpBarycentScalarWeights,
        'minInterpPoints': 3,
      },
      'unstinterp': {
        'method': self.unstinterpBarycentWeights,
        'minInterpPoints': 3,
      },
      'inverseD': {
        'method': self.inverseDistanceWeights,
        'minInterpPoints': 3,
      },
      'barycentricScalar': {
        'method': self.barycentricScalarWeights,
        'minInterpPoints': 4,
      },
      'barycentric': {
        'method': self.barycentricWeights,
        'minInterpPoints': 4,
      },
      'pbarycentric': {
        'method': self.projectedBarycentricWeights,
        'minInterpPoints': 4,
        'distanceMethod': 'cartesian',
      },
    }
    assert weightMethod in list(weightMethods.keys()), \
      'InterpolateCartesian: '+weightMethod+' is not one of the available weightMethods'

    self.weightMethod = weightMethods[weightMethod]['method']

    self.minInterpPoints = weightMethods[weightMethod]['minInterpPoints']

    assert nInterpPoints >= self.minInterpPoints, \
      'InterpolateCartesian: '+weightMethod+' weightMethod requires nInterpPoints > '+str(self.minInterpPoints)

    # init distanceMethod
    self.distanceMethods = {
      'cartesian': self.cartesianInDistances,
    }
    self.distanceMethod = weightMethods[weightMethod].get(
      'distanceMethod', distanceMethod)
    self.distanceMethodInitialized = False

    # init input coordinates
    self.nIn = len(XIn)
    self.XIn = XIn
    self.YIn = YIn
    self.ZIn = ZIn

    XYZ = np.empty((self.nIn, 3))
    for ii in np.arange(0, self.nIn):
      XYZ[ii,:] = np.asarray([self.XIn[ii], self.YIn[ii], self.ZIn[ii]])
    self.Scale = np.mean(np.sqrt(np.square(XYZ).sum(axis=1)))

    self.TreeIn = cKDTree(list(zip(XIn, YIn, ZIn)))
    self.nInterpPoints = nInterpPoints
    self.nnInterpInds = []

  @staticmethod
  def getWorkers(nprocs):
    if nprocs > 1:
      workers = mp.Pool(processes = nprocs)
    else:
      workers = None
    return workers

  @staticmethod
  def finalizeWorkers(workers):
    if workers is not None:
      workers.close()
      workers.join()

  def initNeighbors(self, XOut, YOut, ZOut):
    self.XOut = XOut
    self.YOut = YOut
    self.ZOut = ZOut
    self.nOut = len(XOut)
    self.nnInterpDistances, self.nnInterpInds = \
      self.TreeIn.query(
        list(zip(XOut, YOut, ZOut)),
        k = self.nInterpPoints
      )

  def InDistances(self, Inds1, Inds2):
    if not self.distanceMethodInitialized:
      assert self.distanceMethod in list(self.distanceMethods.keys()), \
        'InterpolateCartesian: '+distanceMethod+' is not one of the available distanceMethods'
      self.InDistanceFunction = self.distanceMethods[self.distanceMethod]
      self.distanceMethodInitialized = True
    return self.InDistanceFunction(Inds1, Inds2)

  def cartesianInDistances(self, Inds1, Inds2):
    # caclculates distances between input nodes for two sets of indices
    # Inds1, Inds2 may be scalars or arrays
    # return cartesian distance
    # + scalar if Inds1 and Inds2 are scalars
    # + otherwise same shape as Inds1 and/or Inds2

    if isinstance(Inds1, Iterable) and isinstance(Inds2, Iterable):
      assert len(Inds1) == len(Inds2), 'Inds1 must be same size as Inds2'

    X1 = self.XIn[Inds1]
    Y1 = self.YIn[Inds1]
    Z1 = self.ZIn[Inds1]

    X2 = self.XIn[Inds2]
    Y2 = self.YIn[Inds2]
    Z2 = self.ZIn[Inds2]

    return CartesianDistance(X1, X2, Y1, Y2, Z1, Z2)

  def initWeights(self, XOut=None, YOut=None, ZOut=None,
              updateNeighbors = False):
    if XOut is not None:
      assert len(XOut) == len(YOut), 'InterpolateCartesian::initWeights: output coordinates must be identical size'
      assert len(XOut) == len(ZOut), 'InterpolateCartesian::initWeights: output coordinates must be identical size'
      if len(self.nnInterpInds) != len(XOut) or updateNeighbors:
        self.initNeighbors(XOut, YOut, ZOut)
    self.weightMethod()

  def inverseDistanceWeights(self):
    self.nnInterpWeights = np.zeros((self.nOut, self.nInterpPoints))
    for kk, nnDistances in enumerate(self.nnInterpDistances):
      if nnDistances.min() < self.eps:
        self.nnInterpWeights[kk,np.argmin(nnDistances)] = 1.0
      else:
        self.nnInterpWeights[kk,:] = 1.0 / nnDistances
        # self.nnInterpWeights[kk,:] = 1.0 / np.square(nnDistances)

    nnInterpWeightSums = np.sum(self.nnInterpWeights, axis=1)
    for k in np.arange(0, self.nnInterpWeights.shape[1]):
      self.nnInterpWeights[:,k] = np.divide(self.nnInterpWeights[:,k], nnInterpWeightSums)

  def unstinterpBarycentScalarWeights(self):
    # calculate weights
    self.nnInterpWeights = np.zeros((self.nOut, self.nInterpPoints))
    bw = np.empty(self.nInterpPoints)
    masks = np.logical_not(np.identity(self.nInterpPoints))
    for kk, (nnInds, nnDistances) in enumerate(
      list(zip(self.nnInterpInds, self.nnInterpDistances))):
      if nnDistances.min()*self.Scale < self.eps:
        self.nnInterpWeights[kk,np.argmin(nnDistances)] = 1.0
      else:
        bw[:] = 0.0
        for jj, (nnInd, nnDistance, mask) in enumerate(
          list(zip(nnInds, nnDistances, masks))):
          otherDistances = self.InDistances(nnInd, nnInds[mask])
          wprod = np.prod(otherDistances)
          bw[jj] = 1.0 / np.prod(otherDistances)

        bsw = 0.
        for jj, (nnInd, nnDistance) in enumerate(
          list(zip(nnInds, nnDistances))):
          bsw += bw[jj] / nnDistance

        self.nnInterpWeights[kk,:] = (bw / nnDistances) / bsw

  def unstinterpBarycentCleanupWeights(self):
    # calculate weights
    self.nnInterpWeights = np.zeros((self.nOut, self.nInterpPoints))
    unstWeights = np.empty(self.nInterpPoints)
    masks = np.logical_not(np.identity(self.nInterpPoints))
    for kk, (nnInds, nnDistances) in enumerate(
      list(zip(self.nnInterpInds, self.nnInterpDistances))):
      if nnDistances.min()*self.Scale < self.eps:
        self.nnInterpWeights[kk,np.argmin(nnDistances)] = 1.0
      else:
        unstWeights[:] = 0.0
        for jj, (nnInd, nnDistance, mask) in enumerate(
          list(zip(nnInds, nnDistances, masks))):
          otherDistances = self.InDistances(nnInd, nnInds[mask])
          unstWeights[jj] = 1.0 / np.prod(otherDistances) / nnDistance
        self.nnInterpWeights[kk,:] = unstWeights / unstWeights.sum()

  def unstinterpBarycentWeights(self):
    # calculate weights
    unstWeights = np.zeros((self.nOut, self.nInterpPoints))
    allOutInds = np.arange(0, self.nOut)

    # handle target locations that are aligned with an input node
    OnNode = self.nnInterpDistances.min(axis=1)*self.Scale < self.eps
    OnNodeInds = allOutInds[OnNode]
    minInds = np.argmin(self.nnInterpDistances[OnNodeInds,:], axis=1)
    for (kk, minInd) in list(zip(OnNodeInds, minInds)):
      unstWeights[kk, minInd] = 1.0

    # handle target locations that are not aligned with any input nodes
    OffNodeInds = allOutInds[np.logical_not(OnNode)]
    otherNodesMasks = np.logical_not(np.identity(self.nInterpPoints))
    otherNodesMaskInds = []
    allNodesInds = np.arange(0, self.nInterpPoints)
    for otherNodesMask in otherNodesMasks:
      otherNodesMaskInds.append(allNodesInds[otherNodesMask])

    unstWeights[OffNodeInds, :] = 1.0 / self.nnInterpDistances[OffNodeInds,:]

    for jj, otherNodesInds in enumerate(otherNodesMaskInds):
      for otherNodeInd in otherNodesInds:
        otherDistances = self.InDistances(
          self.nnInterpInds[OffNodeInds, jj],
          self.nnInterpInds[OffNodeInds, otherNodeInd]
        )
        unstWeights[OffNodeInds, jj] /= otherDistances

    unstWeightsSum = unstWeights[OffNodeInds].sum(axis=1)
    for node in np.arange(0, unstWeights.shape[1]):
      unstWeights[OffNodeInds,node] /= unstWeightsSum

    self.nnInterpWeights = unstWeights

  def barycentricScalarWeights(self):
    nBaryNodes = 3
    minAspectRatio = 0.1

    # setup permutations of local nearest neighbor masks and indices
    subsetIndsPermutations = []
    for inds in itertools.combinations(
      np.arange(0, self.nInterpPoints), nBaryNodes):
      subsetIndsPermutations.append(list(inds))
    subsetIndsPermutations = np.asarray(subsetIndsPermutations)

    # calculate weights
    self.nnInterpWeights = np.zeros((self.nOut, self.nInterpPoints))
    baryWeights = np.empty(self.nInterpPoints)
    for kk, (nnInds, nnDistances) in enumerate(
      list(zip(self.nnInterpInds, self.nnInterpDistances))):

        # loop over subsets of nn nodes until two conditions are met for this target location
        # (1) those nodes form a non-zero area triangle
        # (2) those nodes circumscribe the target location
        for orderedSubsetInds in subsetIndsPermutations:

          # Advance the triangle indices by 1
          # Reduces spurious errors in coarse-to-fine-to-coarse test, by not sure why
          # TODO: check if errors are caused by using great circle instead of planar distances
          subsetInds = np.roll(orderedSubsetInds, 1)

          baryWeights[:] = 0.0
          nnSubsetInds = nnInds[subsetInds]

          # Define a triangle with first two input nodes as the base, third node as the apex
          # + the origin is located at the first node, (x1, y1)
          # + (x1,y1), (x2,y2), (x3,y3) are the coordinates of the vertices
          # + d1, d2, d3 are the lengths of the sides of the triangle
          # + l1, l2, l3 are the distances from the vertices to the target location, (x,y)
          #
          #         n3
          #         /`.
          #        /   `.
          #       /      `. d1
          #   d2 /         `.
          #     /   (x,y)    `.
          #    /      .        `.
          #   /..................`.
          # n1         d3          n2

          d1 = self.InDistances(nnSubsetInds[1], nnSubsetInds[2])
          d2 = self.InDistances(nnSubsetInds[0], nnSubsetInds[2])
          d3 = self.InDistances(nnSubsetInds[0], nnSubsetInds[1])
          x1, y1 = 0., 0.
          x2, y2 = d3, 0.
          x3 = (-d1**2 + d2**2 + d3**2) / (2 * d3)
          radicand = d1**2 - (d3 - x3)**2
          #skip skinny/zero-area triangles with y3 < minAspectRatio * d1
          # possibly a too-conservative restriction for randomly located mesh points,
          # but fine for an unstructured mesh that is carefully ordered
          if (radicand/(d1**2)) < np.square(minAspectRatio):
            continue
          y3 = np.sqrt(radicand)

          # ensure denominator (determinant of transformation matrix) is > 0
          denom = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
          if np.abs(denom) == 0.0: continue

          l1, l2, l3 = nnDistances[subsetInds]

          x = (-l2**2 + l1**2 + d3**2) / (2 * d3)
          radicand = l2**2 - (d3 - x)**2
          if (radicand/(l2**2)) < self.eps:
            # this case indicates simple linear interpolation between the first two nodes
            y = 0.
          else:
            y = np.sqrt(radicand)

          baryWeights[subsetInds[0]] = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / denom
          baryWeights[subsetInds[1]] = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / denom
          baryWeights[subsetInds[2]] = 1. - baryWeights[subsetInds[0]] \
                                          - baryWeights[subsetInds[1]]
          # keep the first triangle that circumscribes the target location,
          # defined as having 0.0 ≤ weights ≤ 1.0
          minBaryCheck = baryWeights.min() >= 0.0
          maxBaryCheck = baryWeights.max() <= 1.0 and baryWeights.max() > 0.0
          if minBaryCheck and maxBaryCheck: break

        assert minBaryCheck and maxBaryCheck, \
'''
  baryWeights checks not satisfied; try increasing nInterpPoints
  baryWeights = '''+str(baryWeights)+'''
  nnInterpDistances = '''+str(nnDistances)+'''
  nnInds = '''+str(nnInds)+'''
  denom = '''+str(denom)+'''
  y = '''+str(y)+'''
  y1 = '''+str(y1)+'''
  y2 = '''+str(y2)+'''
  y3 = '''+str(y3)+'''
  x = '''+str(x)+'''
  x1 = '''+str(x1)+'''
  x2 = '''+str(x2)+'''
  x3 = '''+str(x3)+'''
  l1 = '''+str(l1)+'''
  l2 = '''+str(l2)+'''
  l3 = '''+str(l3)+'''
  d1 = '''+str(d1)+'''
  d1^2 = '''+str(d1**2)+'''
  (d3-x3)^2 = '''+str((d3 - x3)**2)+'''
  l2^2 = '''+str(l2**2)+'''
  (d3-x)^2 = '''+str((d3 - x)**2)+'''
  d2 = '''+str(d2)+'''
  d3 = '''+str(d3)+'''
  LatIn = '''+str(self.LatIn[nnInds]*180.0/np.pi)+'''
  LonIn = '''+str(self.LonIn[nnInds]*180.0/np.pi)+'''
  LatIn = '''+str(self.LatIn[nnSubsetInds]*180.0/np.pi)+'''
  LonIn = '''+str(self.LonIn[nnSubsetInds]*180.0/np.pi)+'''
'''
        self.nnInterpWeights[kk,:] = baryWeights

  def barycentricWeights(self):
    nBaryNodes = 3
    minAspectRatio = 0.01

    # Notes:
    # equilateral triangles; AspectRatio = sqrt(3)/2
    # square grid cells, isosceles triangles; AspectRatio = sqrt(2)/2

    baryWeights = np.zeros((self.nOut, self.nInterpPoints))
    minBaryCheck = np.full(self.nOut, False)
    maxBaryCheck = np.full(self.nOut, False)
    allOutInds = np.arange(0, self.nOut)

    if self.calculateDiagnostics:
      self.combinationUsed = np.zeros(self.nOut, dtype=int)
      self.triangleArea = np.zeros(self.nOut)
      self.determinant = np.zeros(self.nOut)
      self.ycoord = np.zeros(self.nOut)
      self.aspectRatio = np.zeros(self.nOut)

    # loop over subsets of nn nodes until two conditions are met for all target locations
    # (1) those nodes form a non-zero area triangle
    # (2) those nodes circumscribe the target location
    # all calculations are vectorized across KeepSearchingInds
    # the number of remaining target locations decreases as iterations proceed
    for ii, indsSet in enumerate(itertools.combinations(
      np.arange(0, self.nInterpPoints), nBaryNodes)):

      orderedSubsetInds = np.array(list(indsSet))

      # Advance the triangle indices by 1
      # Reduces spurious errors in coarse-to-fine-to-coarse test, by not sure why
      # TODO: check if errors are caused by using great circle instead of planar distances
      subsetInds = np.roll(orderedSubsetInds, 1)

      KeepSearching = np.logical_or(
        np.logical_not(minBaryCheck),
        np.logical_not(maxBaryCheck)
      )
      KeepSearchingInds = allOutInds[KeepSearching]
      nTargets = len(KeepSearchingInds)
      if nTargets == 0: break
      if self.calculateDiagnostics:
        self.combinationUsed[KeepSearchingInds] = ii+1
      #print(nTargets, subsetInds)

      nnSubsetInds = self.nnInterpInds[KeepSearchingInds,:]
      nnSubsetInds = nnSubsetInds[:,subsetInds]

      baryWeights[KeepSearchingInds,:] = 0.0

      # Define a triangle with first two input nodes as the base, third node as the apex
      # + the origin is located at n1 = (x1, y1)
      # + (x1,y1), (x2,y2), (x3,y3) are the coordinates of the vertices
      # + d1, d2, d3 are the lengths of the sides of the triangle
      # + l1, l2, l3 are the distances from the vertices to the target location, (x,y)
      #
      #         n3
      #         /`.
      #        /   `.
      #       /      `. d1
      #   d2 /         `.
      #     /   (x,y)    `.
      #    /      .        `.
      #   /..................`.
      # n1         d3          n2

      d1 = self.InDistances(nnSubsetInds[:,1], nnSubsetInds[:,2])
      d2 = self.InDistances(nnSubsetInds[:,0], nnSubsetInds[:,2])
      d3 = self.InDistances(nnSubsetInds[:,0], nnSubsetInds[:,1])
      x1, y1 = np.zeros(nTargets), np.zeros(nTargets)
      x2, y2 = d3, np.zeros(nTargets)
      x3 = (-d1**2 + d2**2 + d3**2) / (2 * d3)
      radicand = d1**2 - (d3 - x3)**2

      #skip skinny/zero-area triangles with y3 < minAspectRatio * d3
      # possibly a too-conservative restriction for randomly located mesh points,
      # but fine for an unstructured mesh that is carefully ordered
      AspectCheck = (radicand/(d3**2)) >= np.square(minAspectRatio)
      y3 = np.zeros(nTargets)
      y3[AspectCheck] = np.sqrt(radicand[AspectCheck])

      # ensure denominator (determinant of transformation matrix) is > 0
      denom = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
      DeterminantCheck = np.abs(denom) > 0.0

      # only continue calculations for passing triangles
      KeepGoing = np.logical_and(AspectCheck, DeterminantCheck)
      if KeepGoing.sum() == 0: continue

      d3 = d3[KeepGoing]
      x1, y1 = x1[KeepGoing], y1[KeepGoing]
      x2, y2 = x2[KeepGoing], y2[KeepGoing]
      x3, y3 = x3[KeepGoing], y3[KeepGoing]
      denom = denom[KeepGoing]
      KeepSearchingInds = KeepSearchingInds[KeepGoing]
      nTargets = len(KeepSearchingInds)

      if self.calculateDiagnostics:
        self.triangleArea[KeepSearchingInds] = 0.5 * d3 * y3
        self.aspectRatio[KeepSearchingInds] = y3 / d3

      l1 = self.nnInterpDistances[KeepSearchingInds,subsetInds[0]]
      l2 = self.nnInterpDistances[KeepSearchingInds,subsetInds[1]]

      x = (-l2**2 + l1**2 + d3**2) / (2 * d3)
      radicand = l2**2 - (d3 - x)**2

      # y<0 is non-physical
      # y==0 signifies linear interpolation between the first two nodes
      ValidY = (radicand >= 0.0)
      if ValidY.sum() == 0: continue

      y = np.zeros(nTargets)
      y[ValidY] = np.sqrt(radicand[ValidY])

      w1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / denom
      w2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / denom
      w1 = w1[ValidY]
      w2 = w2[ValidY]
      KeepSearchingInds = KeepSearchingInds[ValidY]

      if self.calculateDiagnostics:
        self.determinant[KeepSearchingInds] = denom[ValidY]
        self.ycoord[KeepSearchingInds] = y[ValidY]

      baryWeights[KeepSearchingInds,subsetInds[0]] = w1
      baryWeights[KeepSearchingInds,subsetInds[1]] = w2
      baryWeights[KeepSearchingInds,subsetInds[2]] = 1. - (w1 + w2)

      # keep the first triangle that circumscribes the target location,
      # defined as having 0.0 ≤ weights ≤ 1.0
      for kk in KeepSearchingInds:
        minBaryCheck[kk] = baryWeights[kk,:].min() >= 0.0
        maxBaryCheck[kk] = baryWeights[kk,:].max() <= 1.0 and baryWeights[kk,:].max() > 0.0

    allChecks = np.logical_and(minBaryCheck, maxBaryCheck)

    assert np.all(allChecks), \
'''
  baryWeights checks not satisfied for '''+str(self.nOut - allChecks.sum())+' of '+str(self.nOut)+''' target locations
  --> try increasing nInterpPoints
'''
    self.nnInterpWeights = baryWeights

  def projectedBarycentricWeights(self):
    nBaryNodes = 3
    minAspectRatio = 0.01

    # Notes:
    # equilateral triangles; AspectRatio = sqrt(3)/2
    # square grid cells, isosceles triangles; AspectRatio = sqrt(2)/2

    baryWeights = np.zeros((self.nOut, self.nInterpPoints))
    minBaryCheck = np.full(self.nOut, False)
    maxBaryCheck = np.full(self.nOut, False)
    allOutInds = np.arange(0, self.nOut)

    if self.calculateDiagnostics:
      self.combinationUsed = np.zeros(self.nOut, dtype=int)
      self.triangleArea = np.zeros(self.nOut)
      self.determinant = np.zeros(self.nOut)
      self.ycoord = np.zeros(self.nOut)
      self.aspectRatio = np.zeros(self.nOut)

    # loop over subsets of nn nodes until two conditions are met for all target locations
    # (1) those nodes form a non-zero area triangle
    # (2) those nodes circumscribe the target location
    # all calculations are vectorized across KeepSearchingInds
    # the number of remaining target locations decreases as iterations proceed
    for ii, indsSet in enumerate(itertools.combinations(
      np.arange(0, self.nInterpPoints), nBaryNodes)):

      orderedSubsetInds = np.array(list(indsSet))

      #Note: nearest neighbor index rolling is not needed when the target location is projected
      #      to the plane containing the neighbor nodes (triangle)
      subsetInds = orderedSubsetInds

      KeepSearching = np.logical_or(
        np.logical_not(minBaryCheck),
        np.logical_not(maxBaryCheck)
      )
      KeepSearchingInds = allOutInds[KeepSearching]
      nTargets = len(KeepSearchingInds)
      if nTargets == 0: break
      if self.calculateDiagnostics:
        self.combinationUsed[KeepSearchingInds] = ii+1
      #print(nTargets, subsetInds)

      nnSubsetInds = self.nnInterpInds[KeepSearchingInds,:]
      nnSubsetInds = nnSubsetInds[:,subsetInds]

      baryWeights[KeepSearchingInds,:] = 0.0

      # Define a triangle with first two input nodes as the base, third node as the apex
      # + the origin is located at n1 = (x1, y1)
      # + (x1,y1), (x2,y2), (x3,y3) are the coordinates of the vertices
      # + d1, d2, d3 are the lengths of the sides of the triangle
      # + l1, l2, l3 are the distances from the vertices to the target location, (x,y)
      #
      #         n3
      #         /`.
      #        /   `.
      #       /      `. d1
      #   d2 /         `.
      #     /   (x,y)    `.
      #    /      .        `.
      #   /..................`.
      # n1         d3          n2

      d1 = self.InDistances(nnSubsetInds[:,1], nnSubsetInds[:,2])
      d2 = self.InDistances(nnSubsetInds[:,0], nnSubsetInds[:,2])
      d3 = self.InDistances(nnSubsetInds[:,0], nnSubsetInds[:,1])
      x1, y1 = np.zeros(nTargets), np.zeros(nTargets)
      x2, y2 = d3, np.zeros(nTargets)
      x3 = (-d1**2 + d2**2 + d3**2) / (2 * d3)
      radicand = d1**2 - (d3 - x3)**2

      #skip skinny/zero-area triangles with y3 < minAspectRatio * d3
      # possibly a too-conservative restriction for randomly located mesh points,
      # but fine for an unstructured mesh that is carefully ordered
      AspectCheck = (radicand/(d3**2)) >= np.square(minAspectRatio)
      y3 = np.zeros(nTargets)
      y3[AspectCheck] = np.sqrt(radicand[AspectCheck])

      # ensure denominator (determinant of transformation matrix) is > 0
      denom = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
      DeterminantCheck = np.abs(denom) > 0.0

      # only continue calculations for passing triangles
      KeepGoing = np.logical_and(AspectCheck, DeterminantCheck)
      if KeepGoing.sum() == 0: continue

      d3 = d3[KeepGoing]
      x1, y1 = x1[KeepGoing], y1[KeepGoing]
      x2, y2 = x2[KeepGoing], y2[KeepGoing]
      x3, y3 = x3[KeepGoing], y3[KeepGoing]
      denom = denom[KeepGoing]
      KeepSearchingInds = KeepSearchingInds[KeepGoing]
      nnSubsetInds = nnSubsetInds[KeepGoing,:]
      #nnSubsetInds = nnSubsetInds[np.arange(0, self.nTarget)[KeepGoing],:]
      nTargets = len(KeepSearchingInds)

      if self.calculateDiagnostics:
        self.triangleArea[KeepSearchingInds] = 0.5 * d3 * y3
        self.aspectRatio[KeepSearchingInds] = y3 / d3

      # project XTarget, YTarget, ZTarget onto plane containing
      # n1=(X1, Y1, Z1), n2=(X2, Y2, Z2), and n3=(X3, Y3, Z3)

      XTarget = self.XOut[KeepSearchingInds]
      YTarget = self.YOut[KeepSearchingInds]
      ZTarget = self.ZOut[KeepSearchingInds]

      X1 = self.XIn[nnSubsetInds[:,0]]
      Y1 = self.YIn[nnSubsetInds[:,0]]
      Z1 = self.ZIn[nnSubsetInds[:,0]]

      X2 = self.XIn[nnSubsetInds[:,1]]
      Y2 = self.YIn[nnSubsetInds[:,1]]
      Z2 = self.ZIn[nnSubsetInds[:,1]]

      X3 = self.XIn[nnSubsetInds[:,2]]
      Y3 = self.YIn[nnSubsetInds[:,2]]
      Z3 = self.ZIn[nnSubsetInds[:,2]]

      # define two vectors originating from n1 and extending to n2 and n3
      D3x = X2 - X1
      D3y = Y2 - Y1
      D3z = Z2 - Z1

      D2x = X3 - X1
      D2y = Y3 - Y1
      D2z = Z3 - Z1

      # find a vector normal to the plane, equal to the cross product
      # N = D3 X D2
      Nx, Ny, Nz = np.cross(
        np.array([D3x, D3y, D3z]),
        np.array([D2x, D2y, D2z]),
        axis = 0)

      # normalize the normal vector
      # N = N / ||N||
      NLength = CartesianDistance(Nx, 0, Ny, 0, Nz, 0)
      Nx = Nx / NLength
      Ny = Ny / NLength
      Nz = Nz / NLength

      # determine the distance between XTarget and the plane
      # by projecting (T - n1) onto the plane normal vector, N
      # c = L1 dot N
      L1x = XTarget - X1
      L1y = YTarget - Y1
      L1z = ZTarget - Z1
      c = np.abs(np.array([L1x*Nx, L1y*Ny, L1z*Nz]).sum(axis=0))

      # the projected target vector is equal to the original target vector minus
      # the scaled normal vector from the target vector
      # T' = T - c * N
      XTarget_ = XTarget - c * Nx
      YTarget_ = YTarget - c * Ny
      ZTarget_ = ZTarget - c * Nz

      l1 = CartesianDistance(
        XTarget_, X1,
        YTarget_, Y1,
        ZTarget_, Z1)

      l2 = CartesianDistance(
        XTarget_, X2,
        YTarget_, Y2,
        ZTarget_, Z2)

      x = (-l2**2 + l1**2 + d3**2) / (2 * d3)
      radicand = l2**2 - (d3 - x)**2

      # y<0 is non-physical
      # y==0 signifies linear interpolation between the first two nodes
      ValidY = (radicand >= 0.0)
      if ValidY.sum() == 0: continue

      y = np.zeros(nTargets)
      y[ValidY] = np.sqrt(radicand[ValidY])

      w1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / denom
      w2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / denom
      w1 = w1[ValidY]
      w2 = w2[ValidY]
      KeepSearchingInds = KeepSearchingInds[ValidY]

      if self.calculateDiagnostics:
        self.determinant[KeepSearchingInds] = denom[ValidY]
        self.ycoord[KeepSearchingInds] = y[ValidY]

      baryWeights[KeepSearchingInds,subsetInds[0]] = w1
      baryWeights[KeepSearchingInds,subsetInds[1]] = w2
      baryWeights[KeepSearchingInds,subsetInds[2]] = 1. - (w1 + w2)

      # keep the first triangle that circumscribes the target location,
      # defined as having 0.0 ≤ weights ≤ 1.0
      for kk in KeepSearchingInds:
        minBaryCheck[kk] = baryWeights[kk,:].min() >= 0.0
        maxBaryCheck[kk] = baryWeights[kk,:].max() <= 1.0 and baryWeights[kk,:].max() > 0.0

    allChecks = np.logical_and(minBaryCheck, maxBaryCheck)

    assert np.all(allChecks), \
'''
  baryWeights checks not satisfied for '''+str(self.nOut - allChecks.sum())+' of '+str(self.nOut)+''' target locations
  --> try increasing nInterpPoints
'''
    self.nnInterpWeights = baryWeights

  def apply(self, fIn, OutInds=None):
    assert isinstance(fIn, Iterable) and \
           fIn.shape == self.XIn.shape, 'InterpolateCartesian::apply: fIn must have same shape as original XIn'

    if OutInds is None:
      fOut = np.multiply(
        fIn[self.nnInterpInds],
        self.nnInterpWeights).sum(axis=1)
    elif isinstance(OutInds, Iterable):
      fOut = np.multiply(
        fIn[self.nnInterpInds[OutInds,:]],
        self.nnInterpWeights[OutInds,:]).sum(axis=1)
    else:
      fOut = applyAtIndex(fIn, OutInds)
    return fOut

  def applyAtIndex(self, fIn, OutInd):
    return np.multiply(
      fIn[self.nnInterpInds[OutInd,:]],
      self.nnInterpWeights[OutInd,:]).sum()

class InterpolateLonLat(InterpolateCartesian):
  def __init__(self, LonIn, LatIn, nInterpPoints = 5,
               weightMethod = 'unstinterp',
               distanceMethod = 'greatcircle',
               Radius = 1.,
               calculateDiagnostics = False):
    self.LatIn = LatIn
    self.LonIn = LonIn
    self.Radius = Radius
    XIn, YIn, ZIn = LonLat2Cartesian(LonIn, LatIn, self.Radius)
    super().__init__(XIn, YIn, ZIn, nInterpPoints,
      weightMethod,
      distanceMethod,
      calculateDiagnostics,
    )
    self.Scale = self.Radius
    self.distanceMethods['greatcircle'] = self.greatcircleInDistances

  def initNeighbors(self, LonOut, LatOut, nprocs=1):
    XOut, YOut, ZOut = LonLat2Cartesian(LonOut, LatOut, self.Radius)
    #_logger.info('Initializing Neighbors')
    super().initNeighbors(XOut, YOut, ZOut)

    workers = self.getWorkers(nprocs)

    if self.distanceMethod == 'greatcircle':
      #_logger.info('Calculating great circle distances')
      # find nearest neighbor great circle distances
      if workers is None:
        self.nnInterpDistances = np.empty((len(LatOut), self.nInterpPoints))
        for kk, (nnInds, lon, lat) in enumerate(
          list(zip(self.nnInterpInds, LonOut, LatOut))):
          self.nnInterpDistances[kk, :] = HaversineDistance(
            lon, lat,
            self.LonIn[nnInds], self.LatIn[nnInds],
            self.Radius)
      else:
        distances = []
        for kk, (nnInds, lon, lat) in enumerate(
          list(zip(self.nnInterpInds, LonOut, LatOut))):
          distances.append(workers.apply_async(HaversineDistancePar,
            args = (kk,
              lon, lat,
              self.LonIn[nnInds], self.LatIn[nnInds],
              self.Radius,
            )
          ))

        self.finalizeWorkers(workers)

        self.nnInterpDistances = np.empty((len(LatOut), self.nInterpPoints))
        for d_ in distances:
          d = d_.get()
          self.nnInterpDistances[d['kk'], :] = d['d']


  def initWeights(self, LonOut, LatOut, updateNeighbors = False, nprocs=1):
    if len(self.nnInterpInds) != len(LatOut) or updateNeighbors:
      self.initNeighbors(LonOut, LatOut, nprocs)
    super().initWeights()

  def greatcircleInDistances(self, Inds1, Inds2):
    return HaversineDistance(
      self.LonIn[Inds1], self.LatIn[Inds1],
      self.LonIn[Inds2], self.LatIn[Inds2],
      self.Radius)


def LonLat2Cartesian(lon, lat, R = 1.):
    '''
    calculates cartesion coordinates from longitude and latitude
    of a point on a sphere with radius R

    lon and lat are in radians
    '''
    x =  R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)
    return x,y,z

def CartesianDistance(X1, X2, Y1, Y2, Z1, Z2):
  return np.sqrt((X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2)

def HaversineDistancePar(kk, lon1, lat1, lon2, lat2, R = 1.):
  return {'kk': kk, 'd': HaversineDistance(lon1, lat1, lon2, lat2, R)}

def HaversineDistance(lon1, lat1, lon2, lat2, R = 1.):
    #https://janakiev.com/blog/gps-points-distance-python/
    #lon1, lat1 -- radians, scalar or array
    #lon2, lat2 -- radians, scalar or array
    #R          -- earth radius in meters
    #returns great circle distance(s)
    # + scalar if all inputs are scalars
    # + otherwise same shape as pairs of lat/lon
    if isinstance(lat1, Iterable):
      assert isinstance(lon1, Iterable), 'lon1 must be same type as lat1'
      assert len(lon1) == len(lat1), 'lon1 must be same size as lat1'
    if isinstance(lat2, Iterable):
      assert isinstance(lon2, Iterable), 'lon2 must be same type as lat2'
      assert len(lon2) == len(lat2), 'lon2 must be same size as lat2'
    if isinstance(lat1, Iterable) and isinstance(lat2, Iterable):
      assert len(lat1) == len(lat2), 'lat1 must be same size as lat2'

    dtheta = np.subtract(lat2, lat1)
    dphi   = np.subtract(lon2, lon1)

    a = np.add(np.square(np.sin(np.divide(dtheta, 2.))),
               np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)),
                           np.square(np.sin(np.divide(dphi, 2.)))
               )
        )

    return np.multiply(2.*R, np.arctan2(np.sqrt(a), np.sqrt(1. - a)))
