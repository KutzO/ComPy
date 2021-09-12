from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
# sns.set()
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from matplotlib.gridspec import GridSpec
from PyPDF2 import PdfFileMerger
from .DbManager import DBManager
import shutil
import logging
import multiprocessing
from tqdm import tqdm
import sys
import time 

from .Additional_Software.pyvenn.venn import get_labels, venn2, venn3, venn4, venn5

 
        
###Class for plotting the data previously saved in the database
class PlotTheData():
    def __init__(self, tool, outputpathtmp, outputpath, pathDB, 
                 threads = False, small = 10, medium = 16, large = 27, 
                 dicIDs = False, VcfPlotTable = False, styleyaml = False):
        
        #Define global variables
        self.threads = threads
        self.out = outputpath                               
        self.outtmp = outputpathtmp                         
        self.pathDB = pathDB                                
        self.styleyaml = styleyaml

        plt.style.use(self.styleyaml["plotstyle"])
        self.headersize = self.styleyaml["compare"]["headersize"]
        plt.rcParams[
            "ytick.labelsize"] = self.styleyaml["compare"]["xlabel"]
        plt.rcParams[
            "xtick.labelsize"] = self.styleyaml["compare"]["ylabel"]
        plt.rcParams[
            "axes.labelsize"] = self.styleyaml["compare"]["axlabel"]
        plt.rcParams[
            "axes.titlesize"] = self.styleyaml["compare"]["title"]
        plt.rcParams[
            "legend.fontsize"] = self.styleyaml["compare"]["font"]
        plt.rcParams[
            "legend.title_fontsize"] = self.styleyaml["compare"]["titlefont"]
            
        #Prepare data IDs
        self.dicIDs = dicIDs
        self.dicIDsBam = {}
        self.dicIDsVCF = {}
        if "bam" in self.dicIDs.keys():
            for intID in self.dicIDs["bam"].keys():
                self.dicIDsBam[intID] = self.dicIDs["bam"][intID]
        if "vcf" in self.dicIDs.keys():
            for intID in self.dicIDs["vcf"].keys():
                self.dicIDsVCF[intID] = self.dicIDs["vcf"][intID]
        if len(self.dicIDsVCF.keys()) > 0:
            self.variants = True
        else:
            self.variants = False
   
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        #Initiate the data
        self.GetTheData(VcfPlotTable)
        
        #Start plotting
        pool = multiprocessing.Pool(processes=threads)
        lsJobs = []
        if tool in ["all", "bam", "plot"] and len(self.dicIDsBam.keys()) > 0:
            lsJobs.append(pool.apply_async(self.PlotOnTarget))
            lsJobs.append(pool.apply_async(self.PlotGCvsCov))
            lsJobs.append(pool.apply_async(self.PlotCoverage))
            lsJobs.append(pool.apply_async(self.PlotQCmetrics, args=("1",)))
            lsJobs.append(pool.apply_async(self.PlotQCmetrics, args=("2",)))
 
        for job in tqdm(lsJobs):
            job.get()
        pool.close()
        if tool in ["all", "vcf", "plot"] and self.variants:
            self.PlotVariantData()
        self.compylog.info("Finished plotting")
            
        
        
        
        
        
        
    ###Helper function to extract needed data
    def GetTheData(self, VcfPlotTable):
    
        #Start extracting data from database
        if len(self.dicIDsBam.keys()) > 0:
            
            #Get IDs
            bamcheck = [x for x in self.dicIDsBam.keys()]
            
            #Get names
            self.bamnames = [self.dicIDsBam[x][0] for x in self.dicIDsBam.keys()]
            
            #Get stats
            self.detaildata = DBManager.ExtractData(
                "ReadStatistiks", self.pathDB, ID = bamcheck
            )
            self.detaildata["MeanC"] = pd.to_numeric(self.detaildata["MeanC"])
            self.detaildata["Total"] = pd.to_numeric(self.detaildata["Total"])
            self.detaildata["Mapped"] = pd.to_numeric(self.detaildata["Mapped"])
            self.detaildata["GC"] = pd.to_numeric(self.detaildata["GC"])
            
            #Get QC
            self.phreddata = DBManager.ExtractData(
                "QCmetrics", self.pathDB, ID = bamcheck
            )
            self.phreddata["PHREDmean"] = pd.to_numeric(
                self.phreddata["PHREDmean"]
            )
            self.phreddata["PHREDsd"] = pd.to_numeric(
                self.phreddata["PHREDsd"]
            )
            self.phreddata["ReadCount"] = pd.to_numeric(
                self.phreddata["ReadCount"]
            )
            
            #Get readcount
            self.mappdata = DBManager.ExtractData(
            "ReadMapping", self.pathDB, ID = bamcheck
            )
            self.mappdata["Total"] = pd.to_numeric(self.mappdata["Total"])
            
            #Get info
            self.infobam = DBManager.ExtractData(
                "BamInfo", self.pathDB, ID = bamcheck
            )
            
            #Define Chromosomes and order them to use it in every plotting
            lsChrom = list(set(self.detaildata["Chrom"].values))
            lsGender = []
            if "X" in lsChrom: 
                intIndex = lsChrom.index("X")
                lsGender.append("X")
                lsChrom.pop(intIndex)
            if "Y" in lsChrom:
                intIndex = lsChrom.index("Y")
                lsGender.append("Y")
                lsChrom.pop(intIndex)
            lsChrom = [int(x) for x in lsChrom]
            lsChrom.sort()
            lsChrom += lsGender
            self.lsChrom = [str(x) for x in lsChrom]
        
        #If vcf should also be plotted
        if self.variants:
            vcfcheck = [x for x in self.dicIDsVCF.keys()]
            self.vcfinfo = DBManager.ExtractData(
                "VCFInfo", self.pathDB, ID= vcfcheck    
            )
            self.dfvariants = DBManager.ExtractData(
                "Extracted_Variants", self.pathDB, ID = vcfcheck
            )
            self.dfvariants["Allelefrequency"] = pd.to_numeric(
                self.dfvariants["Allelefrequency"]
            )

        self.vennplots = False
        if VcfPlotTable:
            self.vennplots = pd.read_csv(VcfPlotTable)


    
    ###Helper function to create raw plot grid
    def GetPlotGrid(self, numcols, samples):
        
        #Create main figure
        if len(samples) == 1:
            numrows = 1
            numcols = 1
            intSize = 12
        else:
            numSample = len(samples)
            numrows = numSample // numcols
            if numrows * numcols < numSample:
                numrows += 1
            intSize = 12 * numcols
        imgFigu = plt.figure(
            figsize = (intSize, 6 * numrows),
            constrained_layout = True
        )
        
        #Create outer gridspec
        gridout = GridSpec(
            nrows = numrows, ncols = numcols, figure = imgFigu
        )

        return imgFigu, gridout
    
    
    
    
    
    ###Helper function to plot percentage of reads which are not mapped (per chromosome)
    def PlotOnTarget(self):
        with PdfPages(self.outtmp+"OnTarget.pdf") as pdfFile:
            
            ##Create main figure
            imgFigu, gridout = self.GetPlotGrid(4, self.bamnames)

            plt.suptitle(
                "Percentage of reads on target",
                fontsize = self.headersize
            )
            
            #Start plotting
            cPlt = 0
            for intID in self.dicIDsBam.keys():   
                
                #Get the data
                data = self.detaildata.loc[(self.detaildata["ID"] == intID)]
                datamap = self.mappdata.loc[(self.mappdata["ID"] == intID)]
                datainfo = self.infobam.loc[self.infobam["ID"] == intID]
                
                if datainfo["Reduced"].values[0] == 0:
                    strReduced = "No"
                else:
                    strReduced = "Yes"
                subsamples = datainfo["Subsamples"].values[0]
                fileclass = datainfo["FileClass"].values[0]
                bedid = datainfo["BedID"].values[0]
                flag = datainfo["Flag"].values[0]
                
                lsOnTrg = []
                for chromosome in self.lsChrom:
                    try:
                        lsOnTrg.append(
                            (
                                sum(
                                    data.loc[
                                        data[
                                            "Chrom"
                                            ] == chromosome
                                        ]["Mapped"].to_list()
                            )\
                            /sum(
                                datamap.loc[
                                    datamap[
                                        "Chrom"
                                            ] == chromosome
                                    ]["Total"].to_list()
                                ))\
                            * 100
                        )
                    except Exception as e:
                        lsOnTrg.append(0)
                        self.compylog.info(
                            f"No mapped reads found: {chromosome}"
                        )
                #Plot the prepared data
                subplt = plt.Subplot(imgFigu, gridout[cPlt])
                subplt.bar(
                    x = self.lsChrom, height = lsOnTrg, 
                    color = self.styleyaml["barplots"]["colorcode"]
                )
                subplt.set_ylim(0,100)
                subplt.set_title(
                    str(
                        "% of reads on target \n "
                        +f"Name: {self.dicIDsBam[intID][0]}, "
                        +f"Reduced: {strReduced}, "
                        +f"Subsamples: {subsamples}, "
                        +f"Group: {fileclass}, "
                        +f"BED-ID: {bedid}, "
                        +f"Flag: {flag}" 
                    )
                )
                subplt.set_xlabel("Chromosome")
                subplt.set_ylabel("% reads on target per chromosome")
                imgFigu.add_subplot(subplt)
                cPlt += 1
            
            #Save figure
            pdfFile.savefig(imgFigu)
        
        
        
        
    ###Helper function to plot context of GC content and coverage per chromosome
    def PlotGCvsCov(self):
        with PdfPages(self.outtmp + "LowCov.pdf") as pdfFile:
            
            #Create main figure
            imgFigu, gridout = self.GetPlotGrid(4, self.bamnames)
            plt.suptitle(
                "Targets with coverage below Threshold", 
                fontsize = self.headersize
            )
            
            #Getting data
            lsAdd = []
            lsMax = []
            maxthresh = self.styleyaml["barplots"]["highest"]
            middlethresh = self.styleyaml["barplots"]["middle"]
            lowerthresh = self.styleyaml["barplots"]["lower"]
            lowestthresh = self.styleyaml["barplots"]["lowest"]
            
            for intID in self.dicIDsBam.keys():   
                data = self.detaildata.loc[
                    (self.detaildata["ID"] == intID) \
                    & (self.detaildata["MeanC"] < maxthresh)
                ]
                usedData = pd.DataFrame(
                    columns = (["Chrom", "Thresh", "Number"])
                )
                datainfo = self.infobam.loc[self.infobam["ID"] == intID]
                if datainfo["Reduced"].values[0] == 0:
                    strReduced = "No"
                else:
                    strReduced = "Yes"
                subsamples = datainfo["Subsamples"].values[0]
                fileclass = datainfo["FileClass"].values[0]
                bedid = datainfo["BedID"].values[0]
                flag = datainfo["Flag"].values[0]
                
                for chrom in self.lsChrom:
                    lsAdd.append(
                        [
                            intID, chrom, str(maxthresh), 
                            len(data.loc[(data["Chrom"] == chrom) \
                                         & (data["MeanC"]< maxthresh)].values)
                        ]
                    )
                    lsMax.append(
                        len(data.loc[(data["Chrom"] == chrom) \
                             & (data["MeanC"]< maxthresh)].values)
                    )
                    lsAdd.append(
                        [
                            intID, chrom, str(middlethresh), 
                            len(
                                data.loc[(data["Chrom"] == chrom) \
                                         & (data["MeanC"] < middlethresh)
                                        ].values
                             )
                        ]
                    )
                    lsAdd.append(
                        [
                            intID, chrom, str(lowerthresh), 
                            len(data.loc[(data["Chrom"] == chrom) \
                                         & (data["MeanC"] < lowerthresh)
                                        ].values)
                        ]
                    )
                    lsAdd.append(
                        [
                            intID, chrom, str(lowestthresh), 
                            len(data.loc[(data["Chrom"] == chrom) \
                                         & (data["MeanC"] < lowestthresh)
                                        ].values)
                        ]
                    )
            usedData = pd.DataFrame(
                lsAdd, columns = ["ID", "Chrom", "Thresh", "Number"]
            )
            
            #Start plotting
            cPlt = 0
            for intID in set(usedData["ID"].values):
                ax_sub = imgFigu.add_subplot(gridout[cPlt])
                sns.barplot(
                    y = "Number", x = "Chrom", hue = "Thresh", 
                    palette = self.styleyaml["barplots"]["colorthresh"], 
                    data = usedData.loc[usedData["ID"] == intID], ax = ax_sub
                )
                if self.styleyaml["barplots"]["ymax"] == "max":
                    intMax = max(lsMax)
                else:
                    intMax = self.styleyaml["barplots"]["ymax"]
                ax_sub.set(
                    ylim = (
                        0, intMax
                    )
                )
                ax_sub.set(
                    title="Tragets with mean Coverage below Threshold \n "
                    +f"Name: {self.dicIDsBam[intID][0]}, "
                    +f"Reduced: {strReduced}, "
                    +f"Subsamples: {subsamples}, "
                    +f"Group: {fileclass}, "
                    +f"BED-ID:{bedid}, "
                    +f"Flag: {flag}"
                )     
                ax_sub.set_ylabel("Number of Targets")
                ax_sub.set_xlabel("Chromosome")
                cPlt += 1
            
            #Save figure
            pdfFile.savefig(imgFigu)

        
   

        
    ###Helper function to plot mean per base PHRED score and read length distribution
    def PlotQCmetrics(self, ori):
        with PdfPages(self.outtmp+f"QCphred_{ori}.pdf") as pdfFile:
            
            #Create main figure
            imgFigu, gridout = self.GetPlotGrid(4, self.bamnames)
            plt.suptitle(
                "Mean PHRED score and read length distribution",
                fontsize = self.headersize
            )
            
            #Get the data
            plotdata = self.phreddata.loc[self.phreddata["Mate"] == ori]
            ylimMin = min(self.phreddata["PHREDmean"]) \
                      - max(self.phreddata["PHREDsd"])
            ylimMax = max(self.phreddata["PHREDmean"]) \
                      + max(self.phreddata["PHREDsd"])
            
            #Start plotting
            cPlt = 0
            for intID in self.dicIDsBam.keys():   #One dualplot (PHRED and len) per .bam file
                
                #Get the data
                data = plotdata.loc[plotdata["ID"] == intID]
                datainfo = self.infobam.loc[self.infobam["ID"] == intID]
                if datainfo["Reduced"].values[0] == 0:
                    strReduced = "No"
                else:
                    strReduced = "Yes"
                subsamples = datainfo["Subsamples"].values[0]
                fileclass = datainfo["FileClass"].values[0]
                bedid = datainfo["BedID"].values[0]
                flag = datainfo["Flag"].values[0]
                
                #Normalize the data
                intAll = sum(data["ReadCount"].values)
                lsNorm = [x / intAll for x in data["ReadCount"].values]               
                
                #Plot the data
                ax = imgFigu.add_subplot(gridout[cPlt])
                axes = [ax, ax.twinx()]
                axes[0].plot(
                    data["PHREDmean"].values,
                    color = self.styleyaml["qcplot"]["phredscore"]
                )
                axes[0].fill_between(
                    range(0,len(data["PHREDmean"].values)),
                    (data["PHREDmean"] - data["PHREDsd"]).values,
                    (data["PHREDmean"] + data["PHREDsd"]).values,
                    alpha = 0.3, 
                    color = self.styleyaml["qcplot"]["stderivation"]
                )
                axes[0].set(
                    ylabel = "Mean PHRED with std", 
                    ylim = (ylimMin, ylimMax)
                )
                axes[0].grid(False)
                axes[1].bar(
                    height = lsNorm,
                    x = data["Readposition"].values,
                    color = self.styleyaml["qcplot"]["lengthdist"],
                    alpha = 0.5
                )
                axes[1].set_ylabel(
                    "Read length distribution"    
                )
                axes[1].set_ylim(0, 1)
                axes[1].grid(False)
                ax.set_title(
                    "QC metrics (PHRED and length) \n "
                    +f"Name: {self.dicIDsBam[intID][0]}, "
                    +f"Reduced: {strReduced}, "
                    +f"Subsamples: {subsamples}, "
                    +f"Group: {fileclass}, "
                    +f"BED-ID: {bedid}, "
                    +f"Flag: {flag}"
                )
                ax.legend(
                    ["mean PHRED", "Length distribution"], 
                    labelcolor = ["blue","orange"]
                )
                cPlt += 1
                
            #Save the main plot
            pdfFile.savefig(imgFigu)
                
        
        
        
        
    ###Helper function to plot the mean coverage per target
    def PlotCoverage(self):

        ##Create plot
        with PdfPages(self.outtmp+"CoVplot_Raw.pdf") as pdfFile:
            
            #Create main figure            
            imgFigu, gridout = self.GetPlotGrid(4, self.bamnames)
            plt.suptitle(
                "Mean coverage per targetbase with std", 
                fontsize = self.headersize
            )

            #Start plotting
            cPlt = 0                        
            if self.styleyaml["boxplots"]["ymax"] == "max":
                ymax = max(self.detaildata["MeanC"])
            else:
                ymax = self.styleyaml["boxplots"]["ymax"]
            for intID in self.dicIDsBam.keys():
                
                #Extract from dataframe
                bamdata = self.detaildata.loc[self.detaildata["ID"] == intID]
                datainfo = self.infobam.loc[self.infobam["ID"] == intID]
                if datainfo["Reduced"].values[0] == 0:
                    strReduced = "No"
                else:
                    strReduced = "Yes"
                subsamples = datainfo["Subsamples"].values[0]
                fileclass = datainfo["FileClass"].values[0]
                bedid = datainfo["BedID"].values[0]
                flag = datainfo["Flag"].values[0]
                
                #Plot Coverage
                ax_sub = imgFigu.add_subplot(gridout[cPlt])
                sns.boxplot(
                    y = "MeanC", x = "Chrom", data = bamdata, ax = ax_sub,
                    palette = self.styleyaml["boxplots"]["colorcode"]
                )
                ax_sub.set(
                    title = "MeanC per target on each chromosome \n "
                            +f"Name: {self.dicIDsBam[intID][0]}, "
                            +f"Reduced: {strReduced}, "
                            +f"Subsamples: {subsamples}, "
                            +f"Group: {fileclass}, "
                            +f"BED-ID: {bedid}"
                            +f"Flag: {flag}", 
                   xlabel="Chromosome", 
                   ylabel="Mean coverage"
               )
                ax_sub.grid(False)  
                ax_sub.set_ylim(1, ymax)
                if self.styleyaml["boxplots"]["coverageyscale"] == "log":
                    ax_sub.set_yscale("log")
                cPlt += 1
            
            #Save figure
            pdfFile.savefig(imgFigu)
    
    
  

    
    

            
    ###Main function for variant data plotting
    def PlotVariantData(self):
        
        #Get names and number of files
        self.vcfnames = [self.dicIDsVCF[x][0] for x in self.dicIDsVCF.keys()]
        self.vcfnum = len(self.vcfnames)
        
        #Control if shared variants should be plotted
        if type(self.vennplots) == bool:
            booShare = False
        else:
            booShare = True
        
        # Set minimal samplenumber to 8 to guarantee best plotsize
        # if self.vcfnum < 8:
            # self.vcfnum = 8
        
        #Calculate and plot variant numbers per .vcf file
        lsJobs = []
        pool = multiprocessing.Pool(processes = self.threads)
        for intID in self.dicIDsVCF.keys():
            lsJobs.append(
                pool.apply_async(self.PlotVarAnalysis, args=(intID,))
            )
        for job in tqdm(lsJobs):
            job.get()
        pool.close()
        
        #Plot shared variants 
        if booShare:
            dicComb = {}
            COUNTER = 0
            for comb in self.vennplots.values:
                dicComb[COUNTER] = comb
                COUNTER += 1
            pool = multiprocessing.Pool(self.threads)
            lsJobs = []
            for jobcount in dicComb.keys():
                lsJobs.append(
                    pool.apply_async(
                        self.PlotSharedVariantsInSample,
                        args = (jobcount, dicComb[jobcount])
                    )
                )    
            for job in tqdm(lsJobs):
                job.get()
                
            pool.close()
    
    
    
    
    #Helper function to plot the number of variants found in each .vcf file
    def PlotVarAnalysis(self, intID):
        
        ##Create main figure
        imgFigu, gridout = self.GetPlotGrid(3, [1,2,3]) 
        
        #Set main variables for plotting
        ymax_var_no = max(
            [len(self.dfvariants.loc[self.dfvariants["ID"] == x]) \
             for x in self.dicIDsVCF.keys()]
        )
        intBin = 100 
        if self.styleyaml["barplots"]["vardistyscale"] == "log":
            booLog = True
        else:
            booLog = False
                
        #Start plotting
        cPlt = 0
        
        #Slice dataframe to data of interest
        data = self.dfvariants.loc[self.dfvariants["ID"] == intID]
        
        #Get .vcf info
        datainfo = self.vcfinfo.loc[self.vcfinfo["ID"] == intID]
        fileclass = datainfo["FileClass"].values[0]
        bedid = datainfo["BedID"].values[0]
        
        #Plot AF heat map
        axs0 = imgFigu.add_subplot(gridout[cPlt])
        axs0.set_title(
            str(
                "Number of Variants \n "
                +f"Name: {self.dicIDsVCF[intID][0]}, "
                +f"Group: {fileclass}, "
                +f"BED-ID: {bedid}"
            )
        )                
        sns.histplot(
            x = "Chromosome", y = "Allelefrequency", bins = intBin, 
            data = data, ax = axs0
        )
        cPlt += 1    
        
        #Plot var num
        var_num = len(data)
        ax1 = imgFigu.add_subplot(gridout[cPlt])
        ax1.bar(
            height = var_num, x = f"Variant number : {var_num}",
            color = self.styleyaml["boxplots"]["colorcode"]
        )
        ax1.set(
            title = str(
                "Number of Variants \n "
                +f"Name: {self.dicIDsVCF[intID][0]}, "
                +f"Group: {fileclass}, "
                +f"BED-ID: {bedid}"
            )
        )
        ax1.set_ylim(1, ymax_var_no)
        ax1.set_yscale("log")
        cPlt += 1
        
        #Plot AF dist
        ax2 = imgFigu.add_subplot(gridout[cPlt])
        graph = sns.histplot(
            x = "Allelefrequency",
            data = data,
            bins = intBin,
            alpha = 0.5,
            ax = ax2,
            palette = self.styleyaml["barplots"]["colorcode"],
            log_scale = (False, booLog))
        ax2.set_title(
            str(
                "Allelefrequency distribution of variants \n"
                +f"Name: {self.dicIDsVCF[intID][0]}, "
                +f"Group: {fileclass}, "
                +f"BED-ID: {bedid}"
            )
        )
        cPlt += 1
        pathout = self.outtmp + f"VariantAnalysis_{intID}.pdf"

        #Save figure
        imgFigu.savefig(pathout)
        


                    
                
    
            
            
    ###Helper function to plot VENN diagramms 
    def PlotSharedVariantsInSample(self, count, comb): 
            
            #Define figure
            with PdfPages(
                    self.outtmp+f"SharedVariantsInSample_{count}.pdf"
                    ) as pdfFile:
                
                #Get data
                lsAllVars = []
                lsVCFname = []
                datainfo = self.vcfinfo.loc[self.vcfinfo["ID"].isin(comb)]
                if len(datainfo.values) == 0:
                    self.compylog.error(
                        "No data was found showing comb ID! Maybe wrong IDs?"
                    )
                    raise SystemError(
                        "The IDs in you VCFplot .csv file are not fitting to "
                        + "the input data! Please check!"
                    )
                fileclass = list(set(datainfo["FileClass"]))
                bedid = list(set(datainfo["BedID"]))
                fileclass = datainfo["FileClass"].values[0]
                bedid = datainfo["BedID"].values[0]
                
                #Plot data
                for intID in comb:
                    try:
                        intID == int(intID)
                        data = self.dfvariants.loc[
                            self.dfvariants["ID"] == intID
                        ]
                        plotdata = data.loc[:,
                                            ["Chromosome", "Position", 
                                             "Reference", "Alternative"
                                             ]
                                            ]
                        lsAllVars.append(plotdata.values.tolist())
                        lsVCFname.append(self.dicIDsVCF[intID][0])
                    except:
                        pass
                lsAllVars = [
                    set(tuple([tuple(x) for x in y])) for y in lsAllVars
                ]
                labels = get_labels(lsAllVars)
                if len(lsVCFname) == 2:
                    lsColor = [
                        self.styleyaml["vennplots"]["color1"],
                        self.styleyaml["vennplots"]["color2"]
                    ]
                    fig, ax = venn2(
                        labels, names = lsVCFname, color = lsColor
                    )
                
                elif len(lsVCFname) == 3:
                    lsColor = [
                        self.styleyaml["vennplots"]["color1"],
                        self.styleyaml["vennplots"]["color2"],
                        self.styleyaml["vennplots"]["color3"]
                    ]
                    fig, ax = venn3(
                        labels, names = lsVCFname, color = lsColor
                    )
                
                elif len(lsVCFname) == 4:
                    lsColor = [
                        self.styleyaml["vennplots"]["color1"],
                        self.styleyaml["vennplots"]["color2"],
                        self.styleyaml["vennplots"]["color3"],
                        self.styleyaml["vennplots"]["color4"]
                    ]
                    fig, ax = venn4(
                        labels, names = lsVCFname, color = lsColor
                    )
                
                elif len(lsVCFname) == 5:
                    lsColor = [
                        self.styleyaml["vennplots"]["color1"],
                        self.styleyaml["vennplots"]["color2"],
                        self.styleyaml["vennplots"]["color3"],
                        self.styleyaml["vennplots"]["color4"],
                        self.styleyaml["vennplots"]["color5"]
                    ]
                    fig, ax = venn5(
                        labels, names = lsVCFname, color = lsColor
                    )
                if len(lsVCFname) > 1:
                    ax.set_title(
                        str(
                            "Shared variants of .vcf \n "
                            +f"Group: {fileclass}, "
                            +f"BED-ID: {bedid}"
                        )
                    )
                    
                #Save fig
                pdfFile.savefig(fig)
                plt.close("all")
                    
                    
        
        
        
        
    ###Merge all created plots to create the final output
    def MergeAllPlots(self, dtime):
        
        #Define order of plots
        lsOrder = ["QCphred_1.pdf",
                   "QCphred_2.pdf",
                   "OnTarget.pdf",
                   "CoVplot_Raw.pdf",
                   "LowCov.pdf",
                   "VariantAnalysis",
                   "SharedVariantsInSample"
                   ]
        
        #Collect pathes of all pdfs
        pathes = glob(self.outtmp+"*.pdf")

        #Start collecting plots in order
        pdfFileloc = []
        for filename in lsOrder:
            for i in pathes:
                pathname = i.split("/")[-1]
                if pathname == filename:
                    pdfFileloc.append(i)
            if filename == "VariantAnalysis":
                for intID in self.dicIDsVCF.keys():
                    name = filename + f"_{intID}.pdf"
                    for i in pathes:
                        pathname = i.split("/")[-1]
                        if pathname == name:
                            pdfFileloc.append(i)
            if filename == "SharedVariantsInSample" \
                and type(self.vennplots) != bool:
                    for intcounter in range(len(self.vennplots)):
                        name = filename + f"_{intcounter}.pdf"
                        for i in pathes:
                            pathname = i.split("/")[-1]
                            if name == pathname:
                                pdfFileloc.append(i)
        
        #Merge plots
        Pdfmerge = PdfFileMerger()
        for pdfs in pdfFileloc:
            Pdfmerge.append(pdfs)
        
        #Save merged plot
        Pdfmerge.write(self.out+f"{dtime}_SampleComparison.pdf")
        Pdfmerge.close()
        
        #Delete temp data
        pathdel = self.outtmp[:self.outtmp.find("tmp")+4]
        self.compylog.info(f"Removing tmp folder from {pathdel}")
        shutil.rmtree(pathdel)
        


                
            
        
        