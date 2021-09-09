import pandas as pd  
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from .DbManager import DBManager
import logging
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.table as tbl
import multiprocessing
from tqdm import tqdm
import sys
import statistics
import os


"""
Called by __main__.py if there are more then 3 other samples showing same .bed id.
"""
class PlotCompAll():
    def __init__(self, pathDB, dtime, output, threads, dicIDs, styleyaml):
        
        ###Define global variables
        self.pathDB = pathDB                                                   #Path to current database
        self.dicID = {}
        for k in dicIDs["bam"].keys():
            self.dicID[k] = dicIDs["bam"][k]
        self.newbam = [self.dicID[x][0] for x in self.dicID.keys()]            #New .bam sample names
        self.styleyaml = styleyaml
        plt.style.use(self.styleyaml["plotstyle"])
        #sns.set(style = self.styleyaml["plotstyle"])
        #self.bedid  = bedid                                                    #The .bed identifier of current samples
        #self.subsamples = subsamples                                           #Number of reads collected for QC
        #self.dicID = dicIDs
        #self.reduced = reduced                                                 #True if .bed file was reduced 
        self.out = output
        #if self.reduced:
        #    self.reduced = "Yes"
        #    self.reducedDB = 1
        #else:
        #    self.reduced = "No"
        #    self.reducedDB = 0
        #if FileClass:
        #    self.FileClass = FileClass
        #else:
        #    self.FileClass = "Default"
            
        ###Initialize logging
        self.compylog = logging.getLogger("ComPy")
        
        ###Check if comparison with database is possible (more then 3 other samples with same .bed identifier)
        self.compylog.info("Get Status if comparison is possible")
        
        
        
       
        pool = multiprocessing.Pool(processes=threads)
        lsJobs = []
        ###Start comparison of each new sample with database!
        for intID in self.dicID.keys():
            lsJobs.append(
                pool.apply_async(
                    self.DoPlotting, 
                    args=(
                        intID, dtime,
                    )
                )
            )

        for job in tqdm(lsJobs):
            job.get()
        pool.close()
        self.compylog.info("Finished big comparison of all data")
        
        



        
        
    def DoPlotting(self, intID, dtime):
        
        dfInfo, dfInfoBam, dfInfoVcf, dfInfoBed = self.GetDataInfo(intID)
        dicOldBam = self.CheckSamples(intID, dfInfo) 
        AllStats, AllQC, AllMapp, booComp  = self.GetData(dicOldBam)
        self.compylog.info("Extracted all dataframes from database!")
        if booComp:
            self.compylog.info(
            f"Start plotting big comparison of sample {self.dicID[intID][0]}"
            )
            self.CreatePlotGrid(
                self.dicID[intID][0], len(set(list(AllStats["Chrom"].values)))
            )
            self.compylog.info("Created main frame!")
            self.SubPlotCompare(AllStats, intID, dicOldBam, AllMapp)
            self.compylog.info("Created subplots!")
            self.SubPlotBase(
                intID, AllStats, AllQC, dicOldBam, dfInfoBam, dfInfoVcf, 
                dfInfoBed, AllMapp, dfInfo
            )
            self.SavePlot(intID, dtime)
            self.compylog.info("Saved plots!")
                


    def GetDataInfo(self, intID):
        dfAllInfo = DBManager.ExtractData(
            "BamInfo", self.pathDB, ID = intID
        )
        dfInfoBam = DBManager.ExtractData(
            "BamInfo", self.pathDB, ALL = True
        )
        dfInfoVcf = DBManager.ExtractData(
            "VCFInfo", self.pathDB, ALL = True
        )
        dfInfoBed = DBManager.ExtractData(
            "BedInfo", self.pathDB, ALL = True
        )
        return dfAllInfo, dfInfoBam, dfInfoVcf, dfInfoBed
    
    
    
    
    def CheckSamples(self, intID, dfInfo):
        #Extract all .bam samples from database with current .bed identifier
        bedid = dfInfo["BedID"].values[0]
        if dfInfo["Reduced"].values[0] == 0:
            strReduced = "No"
            strReducedDB = 0
        else:
            strReduced = "Yes"
            strReducedDB = 1
        subsamples = dfInfo["Subsamples"].values[0]
        fileclass = dfInfo["FileClass"].values[0]
        lsAllData = DBManager.ExtractData(
            "BamInfo", self.pathDB, bedid = bedid, 
            strReduce = strReducedDB, subsamples = subsamples,
            FileClass = fileclass
        )        
        #Filter samples if new or in database before
        dicFilteredData = {}
        
        for data in lsAllData:
            booPop = False
            if data[1] == intID:
                booPop = True
            if not booPop:
                dicFilteredData[data[1]] = data[0]
        #Check if there are more then 3 samples left after filtering
        # if len(dicFilteredData.keys()) >= 3:
        #     booCompare = True
        # else:
        #     booCompare = False

        # self.compylog.info(
        #     f"The big comparison is {booCompare} because there are "
        #     +f"{len(dicFilteredData.keys())} data in database to compare"
        # )
        self.compylog.info(
            f"Data from database: {dicFilteredData.values()}"
        )
        # self.compylog.info("New data: {dicFilteredData.values()}")
        self.compylog.info(
            f"Side infos: bedid = {bedid}; "
            +f"subsamples = {subsamples}; "
            +f"reduced= {strReduced}; "
            +f"FileClass = {fileclass}"
        )
        #print(dicFilteredData)
        return dicFilteredData

                 
    def GetData(self, dicOldBam):
        
        #Collect needed information from database
        lsIDs = list(dicOldBam.keys())+list(self.dicID.keys())
        dfStats = DBManager.ExtractData(
            "ReadStatistiks", self.pathDB, ID = lsIDs
        )
        #dfStats["BedID"] = pd.to_numeric(dfStats["BedID"])
        dfStats["Start"] = pd.to_numeric(dfStats["Start"])
        dfStats["End"] = pd.to_numeric(dfStats["End"])
        dfStats["Mapped"] = pd.to_numeric(dfStats["Mapped"])
        dfStats["MeanC"] = pd.to_numeric(dfStats["MeanC"])
        dfStats["GC"] = pd.to_numeric(dfStats["GC"])
        #dfStats["Subsamples"] = pd.to_numeric(dfStats["Subsamples"])
        
        dfQC = DBManager.ExtractData(
            "QCmetrics", self.pathDB, ID = lsIDs
        )
        
        dfMapp = DBManager.ExtractData(
            "ReadMapping", self.pathDB, ID = lsIDs
        )
        dfMapp["Total"] = pd.to_numeric(dfMapp["Total"])
        
        intnumsamples = len(set(dfStats["ID"].values))
        if intnumsamples < 5:
            
            strMsg = str("To few samples to do comparison. Please add more! "
                         + "At least 5 samples are needed, at now only "
                         + f"{intnumsamples} samples are stored in database.."
                     )
            print(strMsg)
            self.compylog.error(strMsg)
            booComp = False
        else:
            booComp = True
        return dfStats, dfQC, dfMapp, booComp
    
    
    
    
    def CreatePlotGrid(self, bamfile, NumChrom):
        imgFigu = plt.figure(figsize = (18, 21))
        imgFigu.suptitle(
            f"Comparison of sample {bamfile} with database \n", size = 20
        )
        self.grid = GridSpec(6, 6, figure = imgFigu)
        self.ax1 = imgFigu.add_subplot(self.grid[0, 0])
        self.ax2 = imgFigu.add_subplot(self.grid[0, 1])
        self.ax3 = imgFigu.add_subplot(self.grid[0, 2])
        self.ax4 = imgFigu.add_subplot(self.grid[1, 0])
        self.ax5 = imgFigu.add_subplot(self.grid[1, 1])
        self.ax6 = imgFigu.add_subplot(self.grid[1, 2])
        self.ax7 = imgFigu.add_subplot(self.grid[0, 3:6])
        self.ax11 = imgFigu.add_subplot(self.grid[1, 3:6])
        self.ax8 = imgFigu.add_subplot(self.grid[2:4, 0:6])
        self.ax9 = imgFigu.add_subplot(self.grid[4:6, 0:3])
        self.ax10 = imgFigu.add_subplot(self.grid[4:6, 3:6])
        plt.subplots_adjust(wspace = 0.7, hspace = 0.7)
        self.imgFigu = imgFigu
        
    
    
    
    def SubPlotBase(self, intID, AllStats, AllQC, dicOldBam, dfInfoBam, 
                    dfInfoVcf, dfInfoBed, AllMapp, dfInfo):
        ###Collect sample metadata
        bedid = dfInfo["BedID"].values[0]
        if dfInfo["Reduced"].values[0] == 0:
            strReduced = "No"
        else:
            strReduced = "Yes"
        subsamples = dfInfo["Subsamples"].values[0]
        fileclass = dfInfo["FileClass"].values[0]
        
        ###Initialize data
        dfUsedDataStats = AllStats.loc[AllStats["ID"] == intID]
        dfUsedDataQC = AllQC.loc[AllQC["ID"] == intID]
        dfUsedDataMap = AllMapp[AllMapp["ID"] == intID]

        ###Bed file information
        numTargets = len(dfUsedDataStats["MeanC"].values)
        meanGC = round(statistics.mean(dfUsedDataStats["GC"]) * 100, 3)
        lsTargetSize = abs(
            dfUsedDataStats["Start"] - dfUsedDataStats["End"]
        ).values
        maxTarget = max(lsTargetSize)
        minTarget = min(lsTargetSize)
        avrgTarget = round(statistics.mean(lsTargetSize), 3)
        
        ###Database information
        numBedfiles = len(dfInfoBed)
        numTotalBam = len(dfInfoBam)
        numTotalVcf = len(dfInfoVcf)
        numGrp = len(set(dfInfoBam["FileClass"]))
        sizeDb = os.path.getsize(self.pathDB)
        booRed = True
        dicSize = {}
        dicSize[0] = "Byte"
        dicSize[1] = "KByte"
        dicSize[2] = "MByte"
        dicSize[3] = "GByte"
        dicSize[4] = "TByte"
        Counter = 0
        while booRed:
            if sizeDb > 1000:
                sizeDb = round(sizeDb / 1000, 3)
                Counter += 1
            else:
                sizeDb = f"{sizeDb} {dicSize[Counter]}"
                booRed = False
        
        ###Calculation information
        numSamples = len(set(AllStats["BAM"].values))
        numUncovTrgt = len([x for x in dfUsedDataStats["MeanC"] if x == 0])
        maxCov = max(dfUsedDataStats["MeanC"].values)
        maxTrgtStart = dfUsedDataStats["Start"].loc[
            dfUsedDataStats["MeanC"] == maxCov
        ].values[0]
        maxTrgtStop = dfUsedDataStats["End"].loc[
            dfUsedDataStats["MeanC"] == maxCov
        ].values[0]
        maxTrgtChrom = dfUsedDataStats["Chrom"].loc[
            dfUsedDataStats["MeanC"] == maxCov
        ].values[0]
        numTotalReads = sum(dfUsedDataMap["Total"].values)
        
        
        
        
        lsTableSampleBedfile = [
            ["ID", intID, "Number targets in bed file", numTargets], 
            ["Bed ID", bedid, "Largest target size", maxTarget], 
            [
                "Reduced targets", strReduced, 
                "Smallest target size", minTarget
            ], 
            [
                "Read samples for QC", subsamples, 
                "Average target size", avrgTarget
            ], 
            ["File group", fileclass, "Average gc content", meanGC], 
        ]
        
        

        lsTableDatabaseCalculation = [
            [
                "Total number of bed files", numBedfiles, 
                "Number of used files", numSamples
            ],
            [
                "Total number of bam", numTotalBam, 
                "Uncovered targets", numUncovTrgt
            ],
            [
                "Total number of vcf", numTotalVcf, 
                "Target with highest cov", 
                f"{maxTrgtChrom}: {maxTrgtStart} - {maxTrgtStop}"
            ],
            [
                "Total number of groups", numGrp, 
                "Total number of reads", numTotalReads
            ],
            [
                "Size database", sizeDb, 
                "-", "-"
            ],
        ]
        table_a = self.ax7.table(
                cellText = lsTableSampleBedfile,
                colLabels = (" ", "Sample info", " ", "Bed file info"),
                bbox = (0, 0, 1, 1),
                # fontsize = (8),
                # colWidths = [0.8]
        )
        self.ax7.axis("off")
        self.ax7.axis("tight")
        self.ax7.set()
        
        strColor = self.styleyaml["tables"]["colorcode"]
        table_a.get_celld()[(0,1)].set_facecolor(strColor)
        table_a.get_celld()[(0,3)].set_facecolor(strColor)
        table_a.get_celld()[(1,0)].set_facecolor(strColor)
        table_a.get_celld()[(2,0)].set_facecolor(strColor)
        table_a.get_celld()[(3,0)].set_facecolor(strColor)
        table_a.get_celld()[(4,0)].set_facecolor(strColor)
        table_a.get_celld()[(5,0)].set_facecolor(strColor)
        table_a.get_celld()[(1,2)].set_facecolor(strColor)
        table_a.get_celld()[(2,2)].set_facecolor(strColor)
        table_a.get_celld()[(3,2)].set_facecolor(strColor)
        table_a.get_celld()[(4,2)].set_facecolor(strColor)
        table_a.get_celld()[(5,2)].set_facecolor(strColor)
        
        table_b = self.ax11.table(
                cellText = lsTableDatabaseCalculation,
                colLabels = (" ", "Database Info", " ", "Calculation Info"),
                bbox = (0, 0, 1, 1),
                # fontsize = (8),
                # colWidths = [0.8]
        )
        table_b.get_celld()[(0,1)].set_facecolor(strColor)
        table_b.get_celld()[(0,3)].set_facecolor(strColor)
        table_b.get_celld()[(1,0)].set_facecolor(strColor)
        table_b.get_celld()[(2,0)].set_facecolor(strColor)
        table_b.get_celld()[(3,0)].set_facecolor(strColor)
        table_b.get_celld()[(4,0)].set_facecolor(strColor)
        table_b.get_celld()[(5,0)].set_facecolor(strColor)
        table_b.get_celld()[(1,2)].set_facecolor(strColor)
        table_b.get_celld()[(2,2)].set_facecolor(strColor)
        table_b.get_celld()[(3,2)].set_facecolor(strColor)
        table_b.get_celld()[(4,2)].set_facecolor(strColor)
        table_b.get_celld()[(5,2)].set_facecolor(strColor)
        
        self.ax11.axis("off")
        self.ax11.axis("tight")
        self.ax11.set()
        # self.ax11.table(
        #     colLabels = ("Bed file info",),
        #     cellText = [
        #         [numSamples], [numTargets], [maxTarget], [minTarget], [avrgTarget],
        #         [meanGC]
        #     ],
        #     rowLabels = [
        #         "Number of files with bed ID", "Number targets in bed file",
        #         "Largest target size", "Smallest target size", 
        #         "Average target size", "Average gc content"
        #     ],
        #     bbox = (1, 0, 1, 1),
        #     fontsize = 8,
        #     colWidths = [0.8]
        # )
        # self.ax11.axis("off")
        # self.ax11.axis("tight")
        # self.ax11.set()
        
        # sns.scatterplot(
            # y="GC", x="MeanC", data=dfUsedDataStats, alpha=0.4, ax=self.ax7
        # )
        # self.ax7.set(
            # title="GC vs mean coverage per .bed file targets"
        # )
        
        ###Plot Coverage
        sns.boxplot(
            y = "MeanC", x = "Chrom", data = dfUsedDataStats, ax = self.ax8, 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        # sns.violinplot(y="MeanC", x="Chrom", data=dfUsedDataStats, ax=self.ax8, cut=0, inner="stick", scale="width")
        # self.ax8.set_ylim(0,max(dfUsedDataStats["MeanC"]+max(dfUsedDataStats["MeanC"]*0.1)))
        self.ax8.set(
            title = "MeanC per target on each chromosome", 
            xlabel = "Chromosome", 
            ylabel = "Mean coverage"
        )
        if self.styleyaml["boxplots"]["coverageyscale"] == "log":
            self.ax8.set(
                yscale = "log"
            )
        self.ax8.grid(False)
        
        ###Plot QC

        #Mate 1
        self.ax9.set(title = "Read QC mate 1", xlabel = "Read nucleotide")
        dfMate1 = dfUsedDataQC.loc[dfUsedDataQC["Mate"] == "1"]
        axes = [self.ax9, self.ax9.twinx()]
        
        axes[0].plot(
            dfMate1["PHREDmean"].values, 
            color = self.styleyaml["qcplot"]["phredscore"]
        )
        axes[0].fill_between(
            range(0, len(dfMate1["PHREDmean"].values)),
            (dfMate1["PHREDmean"] - dfMate1["PHREDsd"]).values,
            (dfMate1["PHREDmean"] + dfMate1["PHREDsd"]).values,
            alpha = 0.3, color = self.styleyaml["qcplot"]["stderivation"]
        )
        axes[0].set(
            ylabel = "Mean PHRED with std",
        )
        axes[0].grid(False)
        intAll = sum(dfMate1["ReadCount"].values)
        lsNorm = [x / intAll for x in dfMate1["ReadCount"].values]
        axes[1].bar(
                    height = lsNorm,
                    x = dfMate1["Readposition"].values,
                    color = self.styleyaml["qcplot"]["lengthdist"],
                    alpha = 0.5
                )
        axes[1].set_ylabel(
            "Read length distribution"    
        )
        axes[1].set_ylim(0, 1)
        axes[1].grid(False)

        
        #Mate 2
        self.ax10.set(title = "Read QC mate 2", xlabel = "Read nucleotide")
        dfMate2 = dfUsedDataQC.loc[dfUsedDataQC["Mate"] == "2"]
        axes = [self.ax10, self.ax10.twinx()]
        axes[0].plot(
            dfMate2["PHREDmean"].values, 
            color = self.styleyaml["qcplot"]["phredscore"]
        )
        axes[0].fill_between(
            range(0, len(dfMate2["PHREDmean"].values)),
            (dfMate2["PHREDmean"] - dfMate2["PHREDsd"]).values,
            (dfMate2["PHREDmean"] + dfMate2["PHREDsd"]).values,
            alpha = 0.3, color = self.styleyaml["qcplot"]["stderivation"]
        )
        axes[0].set(
            ylabel = "Mean PHRED with std",
        )
        axes[0].grid(False)
        intAll = sum(dfMate2["ReadCount"].values)
        lsNorm = [x / intAll for x in dfMate2["ReadCount"].values]
        axes[1].bar(
                    height = lsNorm,
                    x = dfMate2["Readposition"].values,
                    color = self.styleyaml["qcplot"]["lengthdist"],
                    alpha = 0.5
                )
        axes[1].set_ylabel(
            "Read length distribution"    
        )
        axes[1].set_ylim(0, 1)
        axes[1].grid(False)
        
        # axes[1].plot(
        #     dfMate2["ReadCount"].values, 
        #     color = self.styleyaml["qcplot"]["lengthdist"]
        # )
        # axes[1].set(
        #     ylabel = "% of reads showing this length",
        #     ylim = (-0.01, 1)
        # )
        # axes[1].grid(False)
        # self.ax10.legend(
        #     ["mean PHRED", "Length distribution"], 
        #     #labelcolor = ["blue","orange"]
        # )
        
        

            
            
            
            
    def SubPlotCompare(self, dfStats, intID, dicOldBam, dfMapp):
        
        ###Calculate the data
        
        #Slice whole dataframe according to used .bam file 
        dfBamStats = dfStats.loc[dfStats["ID"] == intID]
        dfAllStats = dfStats.loc[dfStats["ID"].isin(dicOldBam.keys())]
        dfMapStatsTarget = dfMapp.loc[dfMapp["ID"] == intID]
        dfMapStatsAll = dfMapp.loc[dfMapp["ID"].isin(dicOldBam.keys())]
        
        #Calculate target size in MB
        dfTargetMB = abs(dfBamStats["Start"] - dfBamStats["End"]) / 1_000_000
        
        #Calculate mean coverage
        intMeanCov = dfBamStats["MeanC"].mean()
        
        #Calculate mean coverage of GC rich regions
        intMeanCovGC = dfBamStats.loc[dfBamStats["GC"] > 0.5]["MeanC"].mean()
        
        #Calculate mean coverage of AT rich regions
        intMeanCovAT = dfBamStats.loc[dfBamStats["GC"] < 0.5]["MeanC"].mean()
    
        #Calculate difference of coverage between GC and AT rich regions
        intDifCov = (intMeanCovGC - intMeanCovAT) / intMeanCov
    
        #Variance of coverage between all targets
        intVarCov = dfBamStats["MeanC"].var()
        lsVarOld = []
        for oldID in dicOldBam.keys():
            lsVarOld.append(
                dfAllStats.loc[dfAllStats["ID"] == oldID]["MeanC"].var()
            )
    
        #Percentage of reads mapped to ROI
        intRoiNew = (
            dfBamStats["Mapped"].sum() / dfMapStatsTarget["Total"].sum()
            )*100
        lsRoiOld = []
        for oldID in dicOldBam.keys():
            lsRoiOld.append(
                (dfStats.loc[
                    dfStats["ID"] == oldID
                        ]["Mapped"].sum() \
                    /dfMapStatsAll.loc[
                        dfMapStatsAll["ID"] == oldID
                        ]["Total"].sum()
                )*100
            )
        
        # #Number of heterozygot positions (2 nucleotide difference and AF >0.1) per MB
        # intNumHet = dfBamStats["Hetero"].sum()/dfTargetMB.sum()
        # lsHetOld = []
        # for oldfile in oldfiles:
            # lsHetOld.append(dfAllStats.loc[dfAllStats["BAM"] == oldfile]["Hetero"].sum()/dfTargetMB.sum())
        

        ###Plot the data
        
        #Mean coverage
        sns.boxplot(
            y = dfAllStats["MeanC"], ax = self.ax1, orient = "vertical", 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        self.ax1.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax1.set_frame_on(False)                    #Makes it transparent
        self.ax1.set(title = "Mean coverage \n", ylabel = "coverage")
        self.ax1.scatter(
            [0], [intMeanCov], s = 200, linewidth = 2, zorder = 10,
            color = self.styleyaml["boxplots"]["markercolor"], 
            marker = self.styleyaml["boxplots"]["markerstyle"],
        )
        self.ax1.grid(False)
        
        #Mean coverage of GC rich regions
        sns.boxplot(
            y = dfAllStats.loc[dfAllStats["GC"] > 0.5]["MeanC"], 
            ax = self.ax2, orient = "vertical", 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        self.ax2.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax2.set_frame_on(False)                    #Makes it transparent
        self.ax2.set(
            title = "Mean coverage \n GC rich regions", ylabel = "coverage"
        )
        self.ax2.scatter(
            [0], [intMeanCovGC], s = 200, linewidth = 2, zorder = 10,
            color = self.styleyaml["boxplots"]["markercolor"], 
            marker = self.styleyaml["boxplots"]["markerstyle"] 
        )
        self.ax2.grid(False)

            
        #Mean coverage of AT rich regions
        sns.boxplot(
            y = dfAllStats.loc[dfAllStats["GC"] < 0.5]["MeanC"], 
            ax = self.ax3, orient = "vertical", 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        self.ax3.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax3.set_frame_on(False)                    #Makes it transparent
        self.ax3.set(
            title = "Mean coverage \n AT rich regions", ylabel = "coverage"
        )
        self.ax3.scatter(
            [0], [intMeanCovAT], s = 200, linewidth = 2, zorder = 10, 
            color = self.styleyaml["boxplots"]["markercolor"], 
            marker = self.styleyaml["boxplots"]["markerstyle"] 
        )
        self.ax3.grid(False)

            
        #Mean coverage variance
        sns.boxplot(
            y = lsVarOld, ax = self.ax4, orient = "vertical", 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        self.ax4.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax4.set_frame_on(False)                    #Makes it transparent
        self.ax4.set(
            title = "Coverage variance \n", ylabel = "variance"
        )
        self.ax4.scatter(
            [0], [intVarCov], s = 200, linewidth = 2, zorder = 10,
            color = self.styleyaml["boxplots"]["markercolor"], 
            marker = self.styleyaml["boxplots"]["markerstyle"] 
        )
        self.ax4.grid(False)
        
        #Mapped to ROI
        sns.boxplot(
            y = lsRoiOld, ax = self.ax5, orient = "vertical", 
            palette = self.styleyaml["boxplots"]["colorcode"]
        )
        self.ax5.axes.get_xaxis().set_visible(False)
        self.ax5.set(title = "Mapped on ROI \n", ylabel = "%")
        self.ax5.scatter(
            [0], [intRoiNew], s = 200, linewidth = 2, zorder = 10,
            color = self.styleyaml["boxplots"]["markercolor"], 
            marker = self.styleyaml["boxplots"]["markerstyle"]
        )
        self.ax5.grid(False)
        
        
        #Coverage pie chart
        tmpfile = self.imgFigu.add_subplot(self.ax6)
        uprthrsh = self.styleyaml["pieplot"]["highest"]
        mdlthrsh = self.styleyaml["pieplot"]["middle"]
        lwrthrsh = self.styleyaml["pieplot"]["lowest"]
        above300 = len(
            [x for x in dfBamStats["MeanC"] if x > uprthrsh]
        )
        above100 = len(
            [x for x in dfBamStats["MeanC"] if x <= uprthrsh and x > mdlthrsh]
        )
        under100 = len(
            [x for x in dfBamStats["MeanC"] if x <= mdlthrsh and x > lwrthrsh]
        )
        under20 = len(
            [x for x in dfBamStats["MeanC"] if x <= lwrthrsh]
        )
        #above50 = len([x for x in dfCov["AN14023706"] if x <= 100 and x > 50])
        #above20 = len([x for x in dfCov["AN14023706"] if x <= 50 and x > 20])
        #below20 = len([x for x in dfCov["AN14023706"] if x <= 20])
    
        totalnum = len([x for x in dfBamStats["MeanC"]])  
        lsPerc = []
        lsPerc.append((above300 / totalnum) * 100)
        lsPerc.append((above100 / totalnum) * 100)
        lsPerc.append((under100 / totalnum) * 100)
        lsPerc.append((under20 / totalnum) * 100)
    
        labels = [
            f"x > {uprthrsh}", 
            f"{uprthrsh} >= x > {mdlthrsh}", 
            f"{mdlthrsh} x > {lwrthrsh}", 
            "x <= {lwrthrsh}"
        ]
        
        tmpfile.pie(
            lsPerc, colors = self.styleyaml["pieplot"]["colorcode"]
        )
        tmpfile.legend(labels, fontsize = 8, bbox_to_anchor = (0.25, 0.25))
        tmpfile.set(title = "Coverage portions")
        

        
    
    
    
    def SavePlot(self, intID, dtime):
            
        #Save plot
        with PdfPages(
                self.out + f"BigCompare/{dtime}_{self.dicID[intID][0]}.pdf"
                ) as pdfFile:
            self.imgFigu.tight_layout()
            pdfFile.savefig(self.imgFigu)
            output = self.out+f"{dtime}_{self.dicID[intID][0]}.pdf"
            self.compylog.info(
                f"Saved plot {output}!"
            )
        
        
        

        
        