#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 10:42:16 2021

@author: kutzolive
"""

###Import needed packages
#External packages
import logging 
import sys

#Own scripts
from .DbManager import DBManager
from .Preparation import DataPreparation
from .PlotData import PlotTheData
from .PlotAll import PlotCompAll

"""
Accepted kwargs: 
    argbamid, argvcfid, argnamelist, argoutput, argdatabase, dtime, argyaml, 
    argThreads
"""
def CompToolPlot(**kwargs):
    ###Prepare logging
    compylog = logging.getLogger("ComPy")
    
    
    """
    1) Get database path
    """
    compylog.info("Prepare database")
    classDataBase = DBManager(kwargs["argDatabase"])    
    

    """
    2) Prepare data
    """
    compylog.info("Prepare arguments")
    classPrep = DataPreparation(
        "plot", lsbamid = kwargs["argbamid"], lsvcfid = kwargs["argvcfid"],
        booDB = classDataBase.booDB, styleyaml = kwargs["argyaml"],
        DBpath= classDataBase.pathDB, nameTable = kwargs["argnamelist"], 
        dtime= kwargs["dtime"], outputpath = kwargs["argoutput"]
    )
    
    """
    3) Plot data
    """
    compylog.info("Plot data")
    print("Create sample comparison")
    PlotClass = PlotTheData(
            "plot", classPrep.outputpathTmp, classPrep.outputpath, 
            classDataBase.pathDB, threads = kwargs["argthreads"], 
            dicIDs = classPrep.dicIDs, styleyaml = classPrep.styleyaml,
            VcfPlotTable = kwargs["argvcfplottable"]
    )
    PlotClass.MergeAllPlots(kwargs["dtime"])

    print("Start individual comparison")
    if kwargs["argbamid"]:
        PlotCompAll(
            classDataBase.pathDB, kwargs["dtime"], classPrep.outputpath, 
            kwargs["argThreads"], classPrep.dicIDs, classPrep.styleyaml,
            fileclasses = kwargs["argclass"]
        )
    else:
        print("No bam files provided.. skipping")
    
    
