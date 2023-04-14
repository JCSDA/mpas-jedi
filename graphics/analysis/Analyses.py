#!/usr/bin/env python3

import inspect
import logging
import multiprocessing as mp

import analysis.StatisticsDatabase as sdb

# AnalysisBase and derived classes
from analysis.AnalysisBase import AnalysisBase

## derived from CategoryBinMethodBase
from analysis.category.CYAxisExpLines import CYAxisExpLines
from analysis.category.FCAxisExpLines import FCAxisExpLines, FCAxisExpLinesDiffCI
from analysis.category.FCScoreCard import FCScoreCard
from analysis.category.CYAxisFCLines import CYAxisFCLines
from analysis.category.CYAxisBinValLines import CYAxisBinValLines

## derived from MultiDimBinMethodBase
from analysis.multidim.CYandBinValAxes2D import CYandBinValAxes2D
from analysis.multidim.FCandBinValAxes2D import FCandBinValAxes2D
from analysis.multidim.BinValAxisProfile import BinValAxisProfile, BinValAxisProfileDiffCI
from analysis.multidim.BinValAxes2D import BinValAxes2D

## derived from AnalysisBase
from analysis.BinValAxisPDF import BinValAxisPDF, BinValAxisPDFMultiExp
from analysis.BinValAxisStatsComposite import BinValAxisStatsComposite
from analysis.GrossValues import GrossValues

AnalysisTypeLookup = {
    'CYAxisExpLines': CYAxisExpLines,
    'FCAxisExpLines': FCAxisExpLines,
    'FCAxisExpLinesDiffCI': FCAxisExpLinesDiffCI,
    'FCScoreCard': FCScoreCard,
    'CYAxisFCLines': CYAxisFCLines,
    'CYAxisBinValLines': CYAxisBinValLines,
    'CYandBinValAxes2D': CYandBinValAxes2D,
    'FCandBinValAxes2D': FCandBinValAxes2D,
    'BinValAxisProfile': BinValAxisProfile,
    'BinValAxisProfileDiffCI': BinValAxisProfileDiffCI,
    'BinValAxes2D': BinValAxes2D,
    'BinValAxisPDF': BinValAxisPDF,
    'BinValAxisPDFMultiExp': BinValAxisPDFMultiExp,
    'BinValAxisStatsComposite': BinValAxisStatsComposite,
    'GrossValues': GrossValues,
}

# NOTES:
# (1) FCAxis* types require non-zero forecast length
# (2) CYAxis* types require > 1 analysis cycle
# (3) CYAxisFCLines requires (1) and (2)
# (4) *DiffCI types require more than one experiment

def AnalysisFactory(db:sdb, analysisType:str, diagnosticGroupings:dict):
    myClass = AnalysisTypeLookup.get(analysisType, None)
    assert (myClass is not None and inspect.isclass(myClass)), (
        '\n\nERROR: AnalysisFactory cannot construct ', analysisType, ' without instructions in AnalysisTypeLookup')
    assert issubclass(myClass, AnalysisBase), (
        '\n\nERROR: AnalysisFactory cannot construct ', analysisType, ', must be subclass of AnalysisBase')
    return myClass(db, analysisType, diagnosticGroupings)


class Analyses():
    '''
    Wrapper class used to initialize and execute derived AnalysisBase classes
    '''
    def __init__(self, db:sdb, analysisTypes:list, diagnosticGroupings:dict = {}, nproc:int = 1):
        self.nproc = nproc
        self.analyses = []
        for anType in analysisTypes:
            self.analyses.append(AnalysisFactory(db, anType, diagnosticGroupings))
        self.logger = logging.getLogger(__name__+'.'+db.DiagSpaceName)
        self.logger.info('Analyses Constructed')

    def analyze(self):
        self.logger.info("Entering Analyses.analyze()")

        if self.nproc > 1:
            workers = mp.Pool(self.nproc)
        else:
            workers = None

        for an in self.analyses:
            an.analyze(workers)

        if workers is not None:
            workers.close()
            workers.join()

        self.logger.info("Exiting Analyses.analyze()")

