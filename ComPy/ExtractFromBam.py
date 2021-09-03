import multiprocessing
import pandas as pd
import pysam
import random
from tqdm import tqdm 
import itertools
import statistics
import math
from .DbManager import DBManager
import logging
import sys

"""
Class for extracting all data from given .bam files

"""
class ExtractInfoData():
    def __init__(self, bam, bamname, ref, threads, dictargets, bedid,
                 pathDB, subsample, version, lsCompatibleVersions, 
                 reducedBed, intID, FileClass, dtime):
        self.version = version                                                 #Current version
        self.lsVersions = lsCompatibleVersions
        self.ID = intID
        self.bam = bam                                                         #.bam file path
        self.ref = ref                                                         #.fasta reference file
        self.threads = threads                      
        self.dicTargets = dictargets                                           #Dictionary with all targets from .bed file (created in script Preparation.py)
        self.bamname = bamname                                                 #Name of the .bam file
        self.bedid = bedid                                                     #The .bed ID (calculated in script Preparation.py)
        self.pathDB = pathDB                                                   #Path to database
        self.subsample = subsample                                             #Number of subsamples taken from every target (default = 100_000)
        self.dtime = dtime
        if reducedBed:                                                         #Define if .bed is reduced (important for data management)
            self.strReduce = 1
        else:
            self.strReduce = 0
        
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
            
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        self.compylog.info("Getting number of reads per chromosome")
        self.GetReadNumbers()
        #self.test = test
            

    ###Extract information about mapped reads, unmapped reads and total number of reads (per chromosome)
    def GetReadNumbers(self):
        ##Building up dataframe with information according to every chromosome in .bed file
        lsDF = []
        with pysam.AlignmentFile(self.bam, "rb") as bamfile:
            for chromosome in self.dicTargets.keys():
                total = [
                    z[3]  for z in bamfile.get_index_statistics() \
                    if z[0] == chromosome
                ][0]
                mapped = [
                    z[1]  for z in bamfile.get_index_statistics() \
                    if z[0] == chromosome
                ][0]
                unmapped = [
                    z[2]  for z in bamfile.get_index_statistics() \
                    if z[0] == chromosome
                ][0]
                lsDF.append(
                    [self.bamname, self.ID, chromosome, total, 
                    mapped, unmapped]
                )
       
        #Transfer the dictionary to a dataframe
        columnnames = [
            "BAM", "ID", "Chrom", "Total", "Mapped", "Unmapped"
        ]
        dfSave = pd.DataFrame(lsDF,columns = columnnames)
        #Adding information about mapped reads to database
        DBManager.InsertData(dfSave.values,"ReadMapping",self.pathDB)

    
    
    
    
    ###Function to determine number of reads taken from every target to reach the number of subsamples per chromosome (default = 100_000)
    def CalcSubsamplesPerTarget(self):
        dicSubsamples = {}
        num_targets = [len(self.dicTargets[x]) for x in self.dicTargets.keys()]
        num_targets = sum(num_targets)
        
        targetsizes = []
        for chromo in self.dicTargets.keys():
            tmpSize = [x[1]-x[0] for x in self.dicTargets[chromo]]
            targetsizes += tmpSize
        
        intTargetTotal = sum(targetsizes)
        trgsizenorm = [
            round((x/intTargetTotal)*self.subsample,0) for x in targetsizes
        ]
        
        if sum(trgsizenorm) > self.subsample:                       
                intDif = sum(trgsizenorm) - self.subsample
                trgsizenorm[trgsizenorm.index(max(trgsizenorm))] -= intDif
        elif sum(trgsizenorm) < self.subsample:
            intDif = self.subsample - sum(trgsizenorm)
            try:
                trgsizenorm[trgsizenorm.index(max(trgsizenorm))] += intDif
            except Exception as e:
                self.compylog.info(
                                "An exception occured while"
                                +" calculating subsamples!"
                                        )
                self.compylog.info(
                                f"Targetsize normalized: {trgsizenorm};"
                                +f" Chromosome: {chromo}; "
                                +f"Number of targets: {num_targets}"
                                        )
                self.compylog.exception(e)
        
        
        for chromo in self.dicTargets.keys():
            dicSubsamples[chromo] = trgsizenorm[:len(self.dicTargets[chromo])]
            trgsizenorm = trgsizenorm[len(self.dicTargets[chromo]):]       

        return dicSubsamples
        
        
        
        
        
    """
    The main extraction function.
    It manages the multiprocessing driven extraction and returns finished lists!
        lsMeanFWD
            Contains mean per base PHRED score of all collected read mate 1
        lsSDEfwd
            Contains standard derivation of mean PHRED saved in lsMeanFWD
        lsMeanREV
            Same as lsMeanFWD but for read mate 2
        lsSDErev
            Same as lsSDEfwd but for read mate 2
        dictLenFwd
            Dictionary containing read length distribution of read mate 1
        dictLenRev
            Same as dictLenFwd but for read mate 2
    Following functions are called:
            GetReadStatistics(), can be found in the same class
            CalcPHRED(), can be found in the same class
    """  
    def ExtractQCdata(self):
    
    
        ##Calculate number of collected reads per target (see function above)
        subsamples = self.CalcSubsamplesPerTarget()
        
        #Define empty output variables (explanation see above this function)
        lsResultPhred = []
        lsResultOri = []
        dictLenFwd = {}
        dictLenRev = {}
        
        #Logging extraction parameters 
        self.compylog.info(f"Reduced .bed file targets: {self.strReduce}")
        self.compylog.info("Number of read samples per chromosome for "
                              +f"QC: {self.subsample}"
        )
        
        
        ##Start extraction  
        pool = multiprocessing.Pool(processes = self.threads)
           
        #Define multiprocessing to extract read statistics + QC data with function GetReadStatistics()
        lsJobs = [pool.apply_async(self.GetReadStatistics, 
                                   args = (
                                       self.dicTargets[chromosom], 
                                       chromosom, 
                                       subsamples[chromosom]
                                   )
                   ) \
                  for chromosom in self.dicTargets.keys()
        ]
        
        
        
        print(f"File: {self.bamname}")
        for job in tqdm(lsJobs):
            RowAdd = []
            chromosom, targets, lsTotal, lsOnTrg, lsmeanc, gc, lsPHRED, \
                lsOrientation, dLenFwd, dLenRev = job.get()
            self.compylog.info(f"Chromosom: {chromosom}")
            for data in zip(targets, lsTotal, lsOnTrg, lsmeanc, gc):
                RowAdd.append(
                                [self.bamname, self.ID, str(chromosom), 
                                 str(data[0][0]), str(data[0][1]), 
                                 str(data[1]), str(data[2]), str(data[3]), 
                                 str(data[4])]
                )
            tmpPHRED = [y[:] for y in lsPHRED]
            tmpOri = [x[:] for x in lsOrientation]
            for data in zip(tmpPHRED, tmpOri):
                lsResultPhred += data[0]
                lsResultOri += data[1]
            for keyname in dLenFwd.keys():
                if keyname in dictLenFwd.keys():
                    dictLenFwd[keyname] += dLenFwd[keyname]
                else:
                    dictLenFwd[keyname] = dLenFwd[keyname]
            for keyname in dLenRev.keys():
                if keyname in dictLenRev.keys():
                    dictLenRev[keyname] += dLenRev[keyname]
                else:
                    dictLenRev[keyname] = dLenRev[keyname]
                    
                    
            #Insert data in table "ReadStatistik"
            DBManager.InsertData(RowAdd,"ReadStatistiks",self.pathDB)
            
            #Add information steppwise to db to avoid large lists/dicts
            del RowAdd
            del tmpPHRED
            del tmpOri
            del lsPHRED
            del lsOrientation
            
        
        
    
        pool.close()
                      #Make sure to lower needed RAM                
        
        
        ##Calculate length distribution and normalize it (max 1 --> 100% of the reads showing the corresponding length)
        # allreadsFwd = sum(dictLenFwd.values())
        # self.compylog.info(
            # "\n Calculating read length distribution from a number of"
            # +f" {allreadsFwd} reads flagged as mate 1!"
        # )
        # print("Calculating read length distribution")
        # for keys in tqdm(dictLenFwd.keys()):
            # dictLenFwd[keys] = dictLenFwd[keys]/allreadsFwd                    #Normalize the data
        # self.compylog.info("Done \n")
        # allreadsRev = sum(dictLenRev.values())
        # self.compylog.info(
            # "\n Calculating read length distribution from a number of "
            # +f"{allreadsRev} reads flagged as mate 2!"
        # )
        # for keys in tqdm(dictLenRev.keys()):
            # dictLenRev[keys] = dictLenRev[keys]/allreadsRev                    #Normalize the data
        # self.compylog.info("Done \n")
        
        
        ##Calculate PHRED scores
        self.compylog.info("Calculating mean PHRED scores!")
        
        #Sort data according to read mate 1 (fwd) and read mate 2 (rev)
        lsPHREDfwd = []
        lsPHREDrev = []
        for _ in zip(lsResultPhred, lsResultOri):
            if _[1] == "1":
                lsPHREDfwd.append(_[0])
            else:
                lsPHREDrev.append(_[0])
        
        #Adding information about read length to define position
        lsPHREDfwd.append([x for x in range(1,max(dictLenFwd.keys())+1)])
        lsPHREDrev.append([x for x in range(1,max(dictLenRev.keys())+1)])
        
        #Merge all PHRED scores from all chromosomes
        zippedphredfwd = list(itertools.zip_longest(*lsPHREDfwd))
        zippedphredrev = list(itertools.zip_longest(*lsPHREDrev))
        
        #Calculate mean + std with CalcPHRED() function
        pool = multiprocessing.Pool(processes=self.threads)
        jobs = [
            pool.apply_async(self.CalcPHRED, args=(Phredrow,)) \
            for Phredrow in zippedphredfwd
        ]                                                                      #only read mate 1
        lsMeanFWD = []
        lsSDEfwd = []
        lsMeanREV = []
        lsSDErev = []
        print("Calculating PHRED scores")
        for _ in tqdm(jobs):
            PhMean, PhSD = _.get()
            if len(PhMean) > 1 and len(PhSD) > 1:
                lsMeanFWD.append(PhMean)
                lsSDEfwd.append(PhSD)
        pool.close()
        
        pool = multiprocessing.Pool(processes=self.threads)
        jobs = [
            pool.apply_async(self.CalcPHRED, args=(Phredrow,)) \
            for Phredrow in zippedphredrev
        ]                                                                      #only read mate 2
        for _ in tqdm(jobs):
            PhMean, PhSD, = _.get()
            if len(PhMean) > 1 and len(PhSD) > 1:
                lsMeanREV.append(PhMean)
                lsSDErev.append(PhSD)
        pool.close()
        
        #Sort the data
        lsMeanFWD.sort(key= lambda x: x[1])
        lsSDEfwd.sort(key = lambda x: x[1])
        lsMeanREV.sort(key= lambda x: x[1])
        lsSDErev.sort(key = lambda x: x[1])

        self.compylog.info("Done! \n")
        return lsMeanFWD, lsSDEfwd, lsMeanREV,lsSDErev, dictLenFwd, dictLenRev





    ###Subprocess to extract needed information from reads
    #Called by ExtractQCdata()
    #Code for GC, threeNucs, hetRet, hetero and meanC made by Dr. rer. nat. Stephan Holger Drukewitz
    #def GetReadStatistics(self, target, chromosom, subsamples, bamfetchtarget):
    def GetReadStatistics(self, targets, chromosom, subsampless):
        
        gc = []
        ##GC content
        # if self.test:
            
        with pysam.FastaFile(self.ref) as fasta:
            for target in targets: 
                seq=fasta.fetch(chromosom,start=target[0]-1, end=target[1]) 
                gc.append(len([x for x in seq if x=='G' or x=='C' ])/len(seq))
        
        ##QC metrics (PHRED)
        lsTargets = []
        lsTotal = []
        lsOnTrg = []
        lsmeanc = []
        arlsReadPhred = []
        arlsReadOri = []
        dicReadLenFwd = {}    
        dicReadLenRev = {}
        with pysam.AlignmentFile(self.bam,"rb") as bamfile:
            for target, subsamples in zip(targets, subsampless):
                lsTargets.append(target)
                #Check total number of reads
                
                
                lsFetchBam = list(
                    bamfile.fetch(
                        chromosom, start=target[0]-1,end=target[1]
                    )
                )
                """
                Iteration durch readklassen und nur mit wichtigen Infos in dataframe
                    --> Aus DF dann random samples ziehen
                """
                
                total_num_reads = len(lsFetchBam)
                
                #Pseudorandom picking of reads
                lsRandom = list(range(0,total_num_reads))
                random.shuffle(lsRandom)
                try:
                    lsRandom = lsRandom[0:int(subsamples)]
                except Exception as e:
                    self.compylog.info(
                        "An error occured at generating random subsamples!"
                    )
                    self.compylog.info(
                        f"Number of subsamples: {int(subsamples)}"
                    )
                    self.compylog.exception(e)
                lsRandom.sort()
                lsRandom.append(0)  #To go on after collecting the last read!
                
                #Start extraction by iterating through all reads 
                
                lsReadPhred = []
                lsReadOri = []
                countreads = 0
                countmatch = 0
                OnTrg = 0
                Total = 0
                for read in lsFetchBam:
                    #Pick rnd reads and determine wether the read is mate 1 ("fwd") or mate 2 ("rev")
                    if len(lsRandom) > 1:
                        if countreads == lsRandom[countmatch]:
                            if read.is_read1:
                                lsReadPhred.append(read.query_qualities)
                                lsReadOri.append("1") 
                                readlength = len(read.query_qualities)
                                if readlength in dicReadLenFwd.keys():
                                    dicReadLenFwd[readlength] += 1
                                else:
                                    dicReadLenFwd[readlength] = 1
                                countmatch += 1
                            elif read.is_read2:
                                lsReadPhred.append(read.query_qualities)
                                lsReadOri.append("2") 
                                readlength = len(read.query_qualities)
                                if readlength in dicReadLenRev.keys():
                                    dicReadLenRev[readlength] += 1
                                else:
                                    dicReadLenRev[readlength] = 1
                                countmatch += 1
                            else:
                                countreads -= 1     #If read is neither flagged as mate 1 nor mate 2 the next read will be picked


                    countreads += 1             #Increase the readcounter
                    
                    #Identify total number of reads and reads which mapped on target
                    Total += 1
                    if read.is_unmapped:
                        pass
                    else:
                        OnTrg += 1
                        
                lsTotal.append(Total)
                lsOnTrg.append(OnTrg)
                lsTargets.append(target)
                arlsReadPhred.append(lsReadPhred)
                arlsReadOri.append(lsReadOri)

                lsLenTotal = []
                for length in dicReadLenFwd.keys():
                    lsLenTotal += [
                        length for x in range(dicReadLenFwd[length])
                    ]
                try:
                    mReadLFwd = sum(lsLenTotal)/len(lsLenTotal)
                except:
                    self.compylog.info(f"No reads found at target: {target}")
                    mReadLFwd = 0
                lsLenTotal = []
                for length in dicReadLenRev.keys():
                    lsLenTotal += [
                        length for x in range(dicReadLenRev[length])
                    ]
                try:
                    mReadLRev = sum(lsLenTotal)/len(lsLenTotal)
                except:
                    mReadLRev = 0
                
                
                meanc = OnTrg *(mReadLFwd + mReadLRev) \
                        / abs(target[0] - target[1])
                lsmeanc.append(meanc)

            return chromosom, targets, lsTotal, lsOnTrg, lsmeanc, gc, \
                    arlsReadPhred, arlsReadOri, dicReadLenFwd, dicReadLenRev

        
        
        
                
    ###Subprocess to calculate mean PHRED and standard derivation
    #Called by ExtractQCdata()
    def CalcPHRED(self, zipline):
        values = [int(x) for x in zipline[:-1] if x != None]
        if len(values) > 1:
            try:
                tmpMean = statistics.mean(values)                                           #Calculate mean
                tmpSDE = math.sqrt(sum(
                                    [(x - tmpMean)**2 for x in values]
                                    ) / (len(values)-1)
                )                                                              #Calculate standard derivation
                #Add result and position in read
                lsMean = [tmpMean,zipline[-1]]                                              
                lsSDE = [tmpSDE,zipline[-1]]
            except Exception as e:
                #print(values)
                print(f"{self.bamname}: tmpMean: {tmpMean}")
                print(e)
                print(
                    f"{self.bamname}: Sum: "
                    +f"{sum([(x - tmpMean)**2 for x in values])}"
                )
                print(f"{self.bamname}: Len: {len(values)-1}")
                sys.exit()
                
        else:
            lsMean = []
            lsSDE = []
        return lsMean, lsSDE





    ###Save QC information and information about mapped reads
    #Called by main program (ComPy.py)
    def SavePhredLen(self, lsMeanPhredFwd, lsStdDevPhredFwd,lsMeanPhredRev, 
                     lsStdDevPhredRev, dicReadLenFwd, dicReadLenRev):   
        
        #Save QC data mate 1
        PosCount = 1
        lsLenValuesFwd = [0 for x in range(0,max(dicReadLenFwd.keys()))]
        for foundlen in dicReadLenFwd.keys():
            lsLenValuesFwd[foundlen-1] = dicReadLenFwd[foundlen]
        data = []
        for entry in zip(lsMeanPhredFwd, lsStdDevPhredFwd, lsLenValuesFwd):
            DBadd = list(entry)
            data.append(
                        [self.bamname, self.ID] \
                        +[PosCount] \
                        +[DBadd[0][0], 
                          DBadd[1][0], DBadd[2], "1"
                          ]
            )
            PosCount += 1
        DBManager.InsertData(data,"QCmetrics",self.pathDB)
        
        #Save QC data mate 2
        PosCount = 1
        lsLenValuesRev = [0 for x in range(0,max(dicReadLenRev.keys()))]
        for foundlen in dicReadLenRev.keys():
            lsLenValuesRev[foundlen - 1] = dicReadLenRev[foundlen]
        data = []
        for entry in zip(lsMeanPhredRev, lsStdDevPhredRev, lsLenValuesRev):
            DBadd = list(entry)
            data.append(
                        [self.bamname, self.ID] \
                        + [PosCount]  \
                        + [DBadd[0][0], 
                           DBadd[1][0], DBadd[2], "2"
                           ]
            )
            PosCount += 1
        DBManager.InsertData(data,"QCmetrics",self.pathDB)
        
        #Adding information that the bam file was processed (and which bed file is used)
        DBManager.SecurityChangeData(
            self.pathDB, self.ID, self.version, "bam"
        )
        

            


        
