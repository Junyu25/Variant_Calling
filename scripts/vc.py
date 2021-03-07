'''
Copyright {2020} Junyu Chen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import os
import argparse
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itertools import repeat
from multiprocessing import Pool, freeze_support


def SimNNK():
    nnkList = []
    protList = []
    for n1 in ["A", "C", "G", "T"]:
        for n2 in ["A", "C", "G", "T"]:
            for k in ["G", "T"]:
                nnk = n1+n2+k
                seq = Seq(nnk)
                mrna = seq.transcribe()
                prot = mrna.translate()
                nnkList.append(nnk)
                protList.append(str(prot))
    return nnkList, protList

def MakeSimDB(fastaFile, OutDir):
    ori = SeqIO.read(fastaFile, "fasta")
    ori = ori.upper()
    mutable_seq = ori.seq.tomutable()
    nnkList, protList = SimNNK()
    seq_list = []
    for i in range(int(len(mutable_seq)/3)):
        mutable_seq = ori.seq.tomutable()
        for nnk, prot in zip(nnkList, protList):
            mutable_seq[i*3:i*3+3] = nnk
            #print(mutable_seq)
            seq = mutable_seq.toseq()
            Seq = SeqRecord(seq)
            Seq.id = str(i+1) + "_" + str(ori.seq[i*3:i*3+3]) + "_" +str(ori.seq[i*3:i*3+3].transcribe().translate()) + "_" + str(i*3+1) + ":" + str(i*3+3) + "_" + nnk + "_" + prot
            Seq.name = ""
            Seq.description = ""
            seq_list.append(Seq)
    outPath = os.path.join(OutDir, "mutate_sim.fasta")
    SeqIO.write(seq_list, outPath, "fasta")
    return outPath

def makeBlastDB(fastaFile, OutDir):
    cmd = "makeblastdb -in " + fastaFile + " -parse_seqids -dbtype nucl -out " + os.path.join(OutDir, os.path.split(fastaFile)[1])
    subprocess.call(cmd, shell=True)
    return os.path.join(OutDir, os.path.split(fastaFile)[1])


def manifestGen(InDir, r1_end, r2_end):
    path = pd.DataFrame()
    for subdir, dirs, files in os.walk(InDir):
        R1 = ""
        R2 = ""
        for file in files:
            if file.endswith(r1_end):
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file.replace(r1_end, r2_end))
                SampleID = file.replace(r1_end, "")
                path = path.append({'SampleID':str(SampleID), "R1":str(R1), "R2":str(R2)}, ignore_index=True)
    return path


## Fastp pair end tirmming and merging
def RunFastp(R1, R2, prefix, OutDir, threads):
    cmd = "fastp --in1 " + R1 + " --in2 " + R2 + \
    " --out1 " + os.path.join(OutDir, prefix + "_R1_unmerged.fastq") + \
    " --out2 " + os.path.join(fastpDir, prefix + "_R2_unmerged.fastq") + \
    " --unpaired1 " + os.path.join(OutDir, prefix + "_R1_unpaired.fastq") + \
    " --unpaired2 " + os.path.join(OutDir, prefix + "_R2_unpaired.fastq") + \
    " --thread " + str(threads) + \
    " --trim_tail1 100 --trim_tail2 100 " + \
    " --merge --merged_out " + os.path.join(OutDir, prefix + ".fastq") + \
    " --correction --overlap_len_require 50" + \
    " --html " + os.path.join(OutDir, prefix + ".html") + \
    " --json " + os.path.join(OutDir, prefix + ".json") + \
    " --report_title " + prefix + "-fastq-merge-report"
    subprocess.call(cmd, shell=True)

## Run fastp in parallel
def RunFastpParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunFastp, zip(R1List, R2List, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## Convert fastq
def fastq2fasta(fastqFile, fastaFile):
    SeqIO.convert(fastqFile, "fastq", fastaFile, "fasta")

def fastq2fastaParallel(fastqFileList, fastaFileList, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(fastq2fasta, zip(fastqFileList, fastaFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunBlastn(fasta, db, OutDir, threads):
    OutFile = os.path.join(OutDir, os.path.split(fasta)[1].replace(".fasta", "") + "_blast.tsv")
    cmd = "blastn -query " + fasta  + " -out " + OutFile + " -evalue 1.0 -max_target_seqs 1 -outfmt 6 -db " + db + " -num_threads " + str(threads) 
    #print(cmd)
    subprocess.call(cmd, shell=True)

def RunBlastnParallel(fastaList, db, OutDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunBlastn, zip(fastaList, repeat(db), repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()



parser = argparse.ArgumentParser(description='Variant Calling')
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-r', '--ref', dest='ref', type=str, required=True,
                    help="the ref path of gene")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='4',
                    help="the number of threads run for a job")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='2',
                    help="the number of jobs run in parallel")
parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='_L001_R1_001.fastq.gz',
                    help="It is the surfix to recognize the forward info, default='_L001_R2_001.fastq.gz'.")
parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='_L001_R2_001.fastq.gz',
                    help="It is the surfix to recognize the reverse info, default='_L001_R2_001.fastq.gz'.")
args = parser.parse_args()

InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
ref = os.path.abspath(args.ref)
threads = int(args.threads)
jobs = int(args.jobs)
r1_end = str(args.sp1)
r2_end = str(args.sp2)

if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)
dbDir = os.path.join(OutDir, "db")
if os.path.exists(dbDir) == 0:
    os.makedirs(dbDir, 0o777, True)
fastpDir = os.path.join(OutDir, "fastp")
if os.path.exists(fastpDir) == 0:
    os.makedirs(fastpDir, 0o777, True)
fastaDir = os.path.join(OutDir, "fasta")
if os.path.exists(fastaDir) == 0:
    os.makedirs(fastaDir, 0o777, True)
blastDir = os.path.join(OutDir, "blast")
if os.path.exists(blastDir) == 0:
    os.makedirs(blastDir, 0o777, True)


path = manifestGen(InDir, r1_end, r2_end)
path.to_csv(os.path.join(OutDir, "manifest.csv"), index=None)
## process manifest
df = path
prefixList = df["SampleID"].tolist()
R1List = df["R1"].tolist()
R2List = df["R2"].tolist()


mutateFile = MakeSimDB(ref, dbDir)
mutateDB = makeBlastDB(mutateFile, dbDir)

#RunFastpParallel(R1List, R2List, prefixList, fastpDir, threads, jobs)
# fastq to fasta
fastqFileList = []
fastaFileList = []
for prefix in prefixList:
    fastqFile = os.path.join(fastpDir, prefix + ".fastq")
    fastaFile = os.path.join(fastaDir, prefix + ".fasta")
    fastqFileList.append(fastqFile)
    fastaFileList.append(fastaFile)

#fastq2fastaParallel(fastqFileList, fastaFileList, jobs)

RunBlastnParallel(fastaFileList, mutateDB, blastDir, threads, jobs)
