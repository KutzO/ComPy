import pandas as pd  
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from .DbManager import DBManager
import logging
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing
from tqdm import tqdm
import sys
"""
Called by __main__.py if there are more then 3 other samples showing same .bed id.
"""
class PlotCompAll():
    def __init__(self, pathDB, newbam, bedid, subsamples, reduced, dtime, 
                 output, threads, dicIDs, FileClass):
        
        ###Define global variables
        self.pathDB = pathDB            #Path to current database
        self.newbam = newbam            #New .bam sample names
        self.bedid  = bedid             #The .bed identifier of current samples
        self.subsamples = subsamples    #Number of reads collected for QC
        self.dicID = dicIDs
        self.reduced = reduced          #True if .bed file was reduced 
        self.out = output
        if self.reduced:
            self.reduced = "Yes"
            self.reducedDB = 1
        else:
            self.reduced = "No"
            self.reducedDB = 0
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
            
        ###Initialize logging
        self.comptoollog = logging.getLogger("ComparisonTool")
        
        ###Check if comparison with database is possible (more then 3 other samples with same .bed identifier)
        self.comptoollog.info("Get Status if comparison is possible")
        booCompare, dicOldBam = self.CheckSamples() 
        
        
       
        pool = multiprocessing.Pool(processes=threads)
        lsJobs = []
        ###Start comparison of each new sample with database!
        for intID in self.dicID.keys():
            lsJobs.append(
                pool.apply_async(
                    self.DoPlotting, 
                    args=(
                        intID, booCompare, dicOldBam, dtime,
                    )
                )
            )

        for job in tqdm(lsJobs):
            job.get()
        pool.close()
        self.comptoollog.info("Finished big comparison of all data")
        
        
        
        
    def CheckSamples(self):
        #Extract all .bam samples from database with current .bed identifier
        lsAllData = DBManager.ExtractData(
            "BamInfo", self.pathDB, bedid = self.bedid, 
            strReduce = self.reducedDB, subsamples = self.subsamples,
            FileClass = self.FileClass
        )        
        #Filter samples if new or in database before
        dicFilteredData = {}
        
        for data in lsAllData:
            booPop = False
            if data[1] in self.dicID.keys():
                booPop = True
            if not booPop:
                dicFilteredData[data[1]] = data[0]
        #Check if there are more then 3 samples left after filtering
        if len(dicFilteredData.keys()) >= 3:
            booCompare = True
        else:
            booCompare = False

        self.comptoollog.info(
            f"The big comparison is {booCompare} because there are "
            +f"{len(dicFilteredData.keys())} data in database to compare"
        )
        self.comptoollog.info(
            f"Data from database: {dicFilteredData.values()}"
        )
        self.comptoollog.info("New data: {dicFilteredData.values()}")
        self.comptoollog.info(
            f"Side infos: bedid = {self.bedid}; "
            +f"subsamples = {self.subsamples}; "
            +f"reduced= {self.reduced}; "
            +f"FileClass = {self.FileClass}"
        )
        
        return booCompare, dicFilteredData



        
        
    def DoPlotting(self, intID, booCompare, dicOldBam, dtime):
        self.comptoollog.info(
            f"Start plotting big comparison of sample {self.dicID[intID][0]}"
        )
        AllStats, AllQC, AllMapp = self.GetData(dicOldBam)
        self.comptoollog.info("Extracted all dataframes from database!")
        self.CreatePlotGrid(
            self.dicID[intID][0], len(set(list(AllStats["Chrom"].values)))
        )
        self.comptoollog.info("Created main frame!")
        if booCompare:
            self.SubPlotCompare(AllStats, intID, dicOldBam, AllMapp)
            self.comptoollog.info("Created subplots!")
        self.SubPlotBase(intID, AllStats, AllQC, dicOldBam)
        self.SavePlot(intID, dtime)
        self.comptoollog.info("Saved plots!")
                


                 
    def GetData(self, dicOldBam):
        
        #Collect needed information from database
        lsIDs = list(dicOldBam.keys())+list(self.dicID.keys())
        dfStats = DBManager.ExtractData(
            "ReadStatistiks", self.pathDB, ID = lsIDs
        )
        dfStats["BedID"] = pd.to_numeric(dfStats["BedID"])
        dfStats["Start"] = pd.to_numeric(dfStats["Start"])
        dfStats["End"] = pd.to_numeric(dfStats["End"])
        dfStats["Mapped"] = pd.to_numeric(dfStats["Mapped"])
        dfStats["MeanC"] = pd.to_numeric(dfStats["MeanC"])
        dfStats["GC"] = pd.to_numeric(dfStats["GC"])
        dfStats["Subsamples"] = pd.to_numeric(dfStats["Subsamples"])
        
        dfQC = DBManager.ExtractData(
            "QCmetrics", self.pathDB, ID = lsIDs
        )
        
        dfMapp = DBManager.ExtractData(
            "ReadMapping", self.pathDB, ID = lsIDs
        )
        dfMapp["Total"] = pd.to_numeric(dfMapp["Total"])
        
        return dfStats, dfQC, dfMapp      
    
    
    
    
    def CreatePlotGrid(self, bamfile, NumChrom):
        imgFigu = plt.figure(figsize=(18,21))
        imgFigu.suptitle(
            f"Comparison of sample {bamfile} with database \n", size=20
        )
        self.grid = GridSpec(6,6, figure=imgFigu)
        self.ax1 = imgFigu.add_subplot(self.grid[0, 0])
        self.ax2 = imgFigu.add_subplot(self.grid[0, 1])
        self.ax3 = imgFigu.add_subplot(self.grid[0, 2])
        self.ax4 = imgFigu.add_subplot(self.grid[1, 0])
        self.ax5 = imgFigu.add_subplot(self.grid[1, 1])
        self.ax6 = imgFigu.add_subplot(self.grid[1, 2])
        self.ax7 = imgFigu.add_subplot(self.grid[0:2, 3:6])
        self.ax8 = imgFigu.add_subplot(self.grid[2:4, 0:6])
        self.ax9 = imgFigu.add_subplot(self.grid[4:6, 0:3])
        self.ax10 = imgFigu.add_subplot(self.grid[4:6, 3:6])
        plt.subplots_adjust(wspace=0.7,hspace=0.7)
        self.imgFigu = imgFigu
        
    
    
    
    def SubPlotBase(self, intID, AllStats, AllQC, dicOldBam):
        
        ###Initialize data
        dfUsedDataStats = AllStats.loc[AllStats["ID"]==intID]
        dfUsedDataQC = AllQC.loc[AllQC["ID"]==intID]
        

        ###Plot GC vs Coverage
        sns.scatterplot(
            y="GC", x="MeanC", data=dfUsedDataStats, alpha=0.4, ax=self.ax7
        )
        self.ax7.set(
            title="GC vs mean coverage per .bed file targets"
        )
        
        ###Plot Coverage
        sns.boxplot(y="MeanC", x="Chrom", data=dfUsedDataStats, ax=self.ax8)
        # sns.violinplot(y="MeanC", x="Chrom", data=dfUsedDataStats, ax=self.ax8, cut=0, inner="stick", scale="width")
        # self.ax8.set_ylim(0,max(dfUsedDataStats["MeanC"]+max(dfUsedDataStats["MeanC"]*0.1)))
        self.ax8.set(
            title="MeanC per target on each chromosome", 
            xlabel="Chromosome", 
            ylabel="Mean coverage"
        )
        self.ax8.grid(False)
        
        ###Plot QC

        #Mate 1
        self.ax9.set(title="Read QC mate 1", xlabel="Read nucleotide")
        dfMate1 = dfUsedDataQC.loc[dfUsedDataQC["Mate"]=="1"]
        axes = [self.ax9, self.ax9.twinx()]
        axes[0].plot(dfMate1["PHREDmean"].values, color="blue")
        axes[0].fill_between(
            range(0,len(dfMate1["PHREDmean"].values)),
            (dfMate1["PHREDmean"]-dfMate1["PHREDsd"]).values,
            (dfMate1["PHREDmean"]+dfMate1["PHREDsd"]).values,
            alpha=0.3
        )
        axes[0].set(ylabel="Mean PHRED with std")
        axes[0].grid(False)
        axes[1].plot(dfMate1["ReadCount"].values, color="orange")
        axes[1].set(ylabel="% of reads showing this length")
        axes[1].grid(False)
        self.ax9.legend(
            ["mean PHRED", "Length distribution"], 
            labelcolor=["blue","orange"]
        )
        
        #Mate 2
        self.ax10.set(title="Read QC mate 2", xlabel="Read nucleotide")
        dfMate2 = dfUsedDataQC.loc[dfUsedDataQC["Mate"]=="2"]
        axes = [self.ax10, self.ax10.twinx()]
        axes[0].plot(dfMate2["PHREDmean"].values, color="blue")
        axes[0].fill_between(
            range(0,len(dfMate2["PHREDmean"].values)),
            (dfMate2["PHREDmean"]-dfMate2["PHREDsd"]).values,
            (dfMate2["PHREDmean"]+dfMate2["PHREDsd"]).values,
            alpha=0.3
        )
        axes[1].plot(dfMate2["ReadCount"].values, color="orange")
        axes[0].set(ylabel="Mean PHRED with std")
        axes[0].grid(False)
        axes[1].set(ylabel="% of reads showing this length")
        axes[1].grid(False)
        self.ax10.legend(
            ["mean PHRED", "Length distribution"], 
            labelcolor=["blue","orange"]
        )
        
        

            
            
            
            
    def SubPlotCompare(self, dfStats, intID, dicOldBam, dfMapp):
        
        ###Calculate the data
        
        #Slice whole dataframe according to used .bam file 
        dfBamStats = dfStats.loc[dfStats["ID"] == intID]
        dfAllStats = dfStats.loc[dfStats["ID"].isin(dicOldBam.keys())]
        dfMapStatsTarget = dfMapp.loc[dfMapp["ID"]==intID]
        dfMapStatsAll = dfMapp.loc[dfMapp["ID"].isin(dicOldBam.keys())]
        
        #Calculate target size in MB
        dfTargetMB = abs(dfBamStats["Start"]-dfBamStats["End"])/1_000_000
        
        #Calculate mean coverage
        intMeanCov = dfBamStats["MeanC"].mean()
        
        #Calculate mean coverage of GC rich regions
        intMeanCovGC = dfBamStats.loc[dfBamStats["GC"]>0.5]["MeanC"].mean()
        
        #Calculate mean coverage of AT rich regions
        intMeanCovAT = dfBamStats.loc[dfBamStats["GC"]<0.5]["MeanC"].mean()
    
        #Calculate difference of coverage between GC and AT rich regions
        intDifCov = (intMeanCovGC-intMeanCovAT)/intMeanCov
    
        #Variance of coverage between all targets
        intVarCov = dfBamStats["MeanC"].var()
        lsVarOld = []
        for oldID in dicOldBam.keys():
            lsVarOld.append(
                dfAllStats.loc[dfAllStats["ID"]==oldID]["MeanC"].var()
            )
    
        #Percentage of reads mapped to ROI
        intRoiNew = (
            dfBamStats["Mapped"].sum()/dfMapStatsTarget["Total"].sum()
            )*100
        lsRoiOld = []
        for oldID in dicOldBam.keys():
            lsRoiOld.append(
                (dfMapStatsAll.loc[
                    dfMapStatsAll["ID"]==oldID
                        ]["Mapped"].sum() \
                    /dfMapStatsAll.loc[
                        dfMapStatsAll["ID"]==oldID
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
        sns.boxplot(y=dfAllStats["MeanC"], ax=self.ax1, orient="vertical")
        self.ax1.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax1.set_frame_on(False)                    #Makes it transparent
        self.ax1.set(title="Mean coverage \n", ylabel="coverage")
        self.ax1.scatter(
            [0], [intMeanCov], color="orange", marker="x", s=200, 
            linewidth=2, zorder=10
        )
        self.ax1.grid(False)
        
        #Mean coverage of GC rich regions
        sns.boxplot(
            y=dfAllStats.loc[dfAllStats["GC"]>0.5]["MeanC"], 
            ax=self.ax2, orient="vertical"
        )
        self.ax2.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax2.set_frame_on(False)                    #Makes it transparent
        self.ax2.set(
            title="Mean coverage \n GC rich regions", ylabel="coverage"
        )
        self.ax2.scatter(
            [0], [intMeanCovGC], color="orange", marker="x", s=200, 
            linewidth=2, zorder=10
        )
        self.ax2.grid(False)

            
        #Mean coverage of AT rich regions
        sns.boxplot(
            y=dfAllStats.loc[dfAllStats["GC"]<0.5]["MeanC"], 
            ax=self.ax3, orient="vertical"
        )
        self.ax3.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax3.set_frame_on(False)                    #Makes it transparent
        self.ax3.set(
            title="Mean coverage \n AT rich regions", ylabel="coverage"
        )
        self.ax3.scatter(
            [0], [intMeanCovAT], color="orange", marker="x", s=200, 
            linewidth=2, zorder=10
        )
        self.ax3.grid(False)

            
        #Mean coverage variance
        sns.boxplot(y=lsVarOld, ax=self.ax4, orient="vertical")
        self.ax4.axes.get_xaxis().set_visible(False)    #Keine x-Axe
        #self.ax4.set_frame_on(False)                    #Makes it transparent
        self.ax4.set(title="Coverage variance \n", ylabel="variance")
        self.ax4.scatter(
            [0], [intVarCov], color="orange", marker="x", s=200, 
            linewidth=2, zorder=10
        )
        self.ax4.grid(False)
        
        #Mapped to ROI
        sns.boxplot(y=lsRoiOld, ax=self.ax5, orient="vertical")
        self.ax5.axes.get_xaxis().set_visible(False)
        self.ax5.set(title="Mapped on ROI \n", ylabel="%")
        self.ax5.scatter(
            [0], [intRoiNew], color="orange", marker = "x", s=200, 
            linewidth=2, zorder=10
        )
        self.ax5.grid(False)
        
        
        #Number of heterozygous positions
        # sns.boxplot(y=lsHetOld, ax = self.ax6, orient= "vertical")
        # self.ax6.axes.get_xaxis().set_visible(False)
        # #self.ax5.set_frame_on(False)
        # self.ax6.set(title="Heterozygous positions \n per MB target", ylabel="Heterozygous positions")
        # self.ax6.scatter([0], [intNumHet], color="orange", marker="x", s=200, linewidth=2, zorder = 10)
        # self.ax6.grid(False)
        
    
    
    
    def SavePlot(self, intID, dtime):
            
        #Save plot
        with PdfPages(
                self.out+f"BigCompare/{dtime}_{self.dicID[intID][0]}.pdf"
                ) as pdfFile:
            self.imgFigu.tight_layout()
            pdfFile.savefig(self.imgFigu)
            output = self.out+f"{dtime}_{self.dicID[intID][0]}.pdf"
            self.comptoollog.info(
                f"Saved plot {output}!"
            )
        
        
        

        
        