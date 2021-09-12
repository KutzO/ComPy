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
import numpy as np

"""
Class for extracting all data from given .bam files
"""
class ExtractInfoData():
    def __init__(self, bam, bamname, ref, threads, dictargets, bedid,
                 pathDB, subsample, version, lsCompatibleVersions, 
                 reducedBed, intID, FileClass, dtime, flag):
        
        ###Define global variables
        self.version = version                                                 
        self.lsVersions = lsCompatibleVersions
        self.ID = intID
        self.bam = bam                                                         
        self.ref = ref                                                         
        self.threads = threads                      
        self.dicTargets = dictargets                                           
        self.bamname = bamname                                                 
        self.bedid = bedid                                                     
        self.pathDB = pathDB                                                   
        self.subsample = subsample                                             
        self.dtime = dtime
        self.flag = flag
        if reducedBed:                                                         
            self.strReduce = 1
        else:
            self.strReduce = 0
        
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
        
        ###Prepare FLAG infos
        self.dicHexa = {
            1   :   1,
            2   :   2,
            4   :   4,
            8   :   8,
            16  :  10,
            32  :  20,
            64  :  40,
            128 :  80,
            256 : 100,
            512 : 200,
            1024: 400,
            2048: 800,
        }

        self.dicFlagFilter = {
            1  : self.paired,
            2  : self.proper,
            4  : self.unmapp,
            8  : self.mateunmap,
            10 : self.revori,
            20 : self.materevori,
            40 : self.readone,
            80 : self.readtwo,
            100: self.secundary,
            200: self.qcfail,
            400: self.dupli,
            800: self.suply,
        }
        
        ###Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        ###Get read numbers
        self.GetReadNumbers()
        
        ###Get used Flag filters
        self.FlagFilter = self.GetFlags(flag)
            
        

    ###Extract information about mapped reads, unmapped reads and total number of reads (per chromosome)
    def GetReadNumbers(self):
        self.compylog.info("Getting number of reads per chromosome")
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
       
        #Transfer the dictionary to a dataframe and save to table "ReadMapping"
        self.compylog.info("Saving read numbers to table ReadMapping")
        columnnames = [
            "BAM", "ID", "Chrom", "Total", "Mapped", "Unmapped"
        ]
        dfSave = pd.DataFrame(lsDF,columns = columnnames)
        DBManager.InsertData(dfSave.values,"ReadMapping",self.pathDB)
        self.compylog.info("Done extracting read numbers!")

    
    
    
    
    ###Function to determine number of reads taken from every target to reach the number of subsamples per chromosome
    def CalcSubsamplesPerTarget(self):
        self.compylog.info(f"Calculate subsamples ({self.subsample} per target")
        self.compylog.info(f"Reduced .bed file targets: {self.strReduce}")
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
        self.compylog.info("Done calculating subsamples per target")
        return dicSubsamples
        
        
    ###Function to translate given flag to flag filter functions
    def GetFlags(self, number):
        self.compylog.info("Define flag filter functions")
        lsFunc = []
        while True:
            lsHexa = list(self.dicHexa.keys())
            for i in range(0, len(lsHexa)):
                value = lsHexa[i]
                if value > number:
                    number -= lsHexa[i-1]
                    lsFunc.append((self.dicHexa[lsHexa[i-1]]))
                    break
            if number > max(lsHexa):
                number -= max(lsHexa)
                lsFunc.append(self.dicHexa[max(lsHexa)])
            if number == 0:
                break
        lsFunc = list(set(lsFunc))
        lsFunc.sort()
        self.compylog.info("Done filtering flags!")
        self.compylog.info(
            f"Flag: {self.flag} means filter {lsFunc}"
        )
        return lsFunc

    """
    Helper functions to filter flags
    """
    def unmapp(self, read):
        if read.is_unmapped:
            return True
        else:
            return False

    def paired(self, read):
        if read.is_paired:
            return True
        else:
            return False

    def proper(self, read):
        if read.is_proper_pair:
            return True
        else:
            return False

    def mateunmap(self, read):
        if read.mate_is_unmapped:
            return True
        else:
            return False

    def revori(self, read):
        if read.is_reverse:
            return True
        else:
            return False

    def materevori(self, read):
        if read.mate_is_reverse:
            return True
        else:
            return False
        
    def readone(self, read):
        if read.is_read1:
            return True
        else:
            return False

    def readtwo(self, read):
        if read.is_read2:
            return True
        else:
            return False

    def secundary(self, read):
        if read.is_secondary:
            return True
        else:
            return False

    def qcfail(self, read):
        if read.is_qcfail:
            return True
        else:
            return False

    def dupli(self, read):
        if read.is_duplicate:
            return True
        else:
            return False

    def suply(self, read):
        if read.is_supplementary:
            return True
        else:
            return False    
    
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
    """  
    def ExtractQCdata(self):
        self.compylog.info("Extract QC data")
    
        ##Calculate number of collected reads per target (see function above)
        subsamples = self.CalcSubsamplesPerTarget()
        
        #Define empty output variables
        lsResultPhred = []
        lsResultOri = []
        dictLenFwd = {}
        dictLenRev = {}
        
        ##Start extraction  
        pool = multiprocessing.Pool(processes = self.threads)
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
        self.compylog.info(f"File: {self.bamname}")
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
            
            #Create big list with ALL extracted QC scores
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
            
            #Delete variables to save storage
            del RowAdd
            del tmpPHRED
            del tmpOri
            del lsPHRED
            del lsOrientation
        pool.close()   
        
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
        lsPHREDfwd.append([x for x in range(1, max(dictLenFwd.keys()) + 1)])
        lsPHREDrev.append([x for x in range(1, max(dictLenRev.keys()) + 1)])
        
        #Merge all PHRED scores from all chromosomes
        zippedphredfwd = list(itertools.zip_longest( * lsPHREDfwd))
        zippedphredrev = list(itertools.zip_longest( * lsPHREDrev))
        
        #Calculate mean + std with CalcPHRED() function
        pool = multiprocessing.Pool(processes=self.threads)
        jobs = [
            pool.apply_async(self.CalcPHRED, args=(Phredrow,)) \
            for Phredrow in zippedphredfwd
        ]                                                                     
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
        ]                                                                      
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

        self.compylog.info("Done Extracting data!")
        return lsMeanFWD, lsSDEfwd, lsMeanREV,lsSDErev, dictLenFwd, dictLenRev





    ###Subprocess to extract needed information from reads
    def GetReadStatistics(self, targets, chromosom, subsampless):

        ##GC content
        gc = []
        with pysam.FastaFile(self.ref) as fasta:
            for target in targets: 
                seq=fasta.fetch(chromosom,start=target[0]-1, end=target[1]) 
                gc.append(len([x for x in seq if x=='G' or x=='C' ])/len(seq))
        
        ##QC metrics (PHRED / LEN) & Coverage
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
                lsFetchBam = list(
                    bamfile.fetch(
                        chromosom, start=target[0]-1,end=target[1]
                    )
                )
                
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
                lsCovLen = []
                dicNames = {}
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
                    
                    #Identify total number of reads and reads passed flag filter and get coverage data
                    Total += 1
                    lsFlagCheck = []
                    for i in self.FlagFilter:
                        lsFlagCheck.append(self.dicFlagFilter[i](read))
                    if True in lsFlagCheck:
                        pass
                    else:
                        OnTrg += 1
                        if read.qname not in dicNames.keys():
                            dicNames[read.qname] = [read.query_alignment_start, read.query_alignment_end]
                        else:
                            if read.query_alignment_start in range(dicNames[read.qname][0], dicNames[read.qname][1]):
                                newend = read.query_alignment_length - (dicNames[read.qname][1] - read.query_alignment_start)
                                realend = dicNames[read.qname][1] + newend
                                dicNames[read.qname][1] = realend
                            else:
                                dicNames[read.qname][1] = dicNames[read.qname][1] + read.query_alignment_length
                
                #Collect data
                lsTotal.append(Total)
                lsOnTrg.append(OnTrg)
                lsTargets.append(target)
                arlsReadPhred.append(lsReadPhred)
                arlsReadOri.append(lsReadOri)
                
                #Calculate coverage
                if OnTrg != 0:
                    meanc = sum([x[1] - x[0] for x in dicNames.values()]) / abs(target[0] - target[1])
                else:
                    self.compylog.info(
                        f"No reads found! {chromosom}, {target}"
                    )
                    meanc = 0
                lsmeanc.append(meanc)

            return chromosom, targets, lsTotal, lsOnTrg, lsmeanc, gc, \
                    arlsReadPhred, arlsReadOri, dicReadLenFwd, dicReadLenRev

        
        
        
                
    ###Subprocess to calculate mean PHRED and standard derivation
    def CalcPHRED(self, zipline):
        
        #Extract QC results without read position
        values = [int(x) for x in zipline[:-1] if x != None]
        if len(values) > 1:
            try:
                #Calculate mean and standard derivation of PHRED scores
                tmpMean = statistics.mean(values)                                           
                tmpSDE = math.sqrt(sum(
                                    [(x - tmpMean)**2 for x in values]
                                    ) / (len(values)-1)
                )            
                
                #Add result and position in read
                lsMean = [tmpMean,zipline[-1]]                                              
                lsSDE = [tmpSDE,zipline[-1]]
            except Exception as e:
                self.compylog.error(f"{self.bamname}: tmpMean: {tmpMean}")
                self.compylog.exception(e)
                self.compylog.error(
                    f"{self.bamname}: Sum: "
                    +f"{sum([(x - tmpMean)**2 for x in values])}"
                )
                self.compylog.error(f"{self.bamname}: Len: {len(values)-1}")
                sys.exit()
                
        else:
            lsMean = []
            lsSDE = []
        return lsMean, lsSDE





    ###Save QC information and information about mapped reads
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
        

            


        
