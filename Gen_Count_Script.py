#!/usr/bin/env python3
import sys, os
import math

def ensure_dir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)

def query_timestamp(site, runNo):
    infodir = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_IBD_0p7MeV/RunInfo/"
    infofile = infodir+f"{site}/{runNo}RunInfo.txt"
    with open(infofile, "r") as f:
        timestamp = int(f.readline().split()[2])
        return timestamp

timestamp_cutingedge = 1388505600

if __name__ == "__main__":

    site = sys.argv[1]

    StreamListDir = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_IBD_0p7MeV/Output/"+site

    basedir = "./counting_result"

    print(basedir)

    cshfiles = ""

    scriptDir = basedir + "Script/" + site + "/"
    outputDir = basedir + "Output/" + site + "/"
    infoDir = basedir + "RunInfo/" + site + "/"
    ensure_dir(scriptDir)
    ensure_dir(outputDir)
    ensure_dir(infoDir)

    ListOfList = os.listdir(StreamListDir)
    for runIdx, run in enumerate(sorted(ListOfList)):
        ext = run[-5:]
        if ext == ".root":

            base = run[:-5]
            runN = run[:5]

            if query_timestamp(site, runN)> timestamp_cutingedge:
                break

            outputfilename = outputDir + runN + ".count.root"
            outputinfoname = infoDir + runN + "RunInfo.txt"
            inputfilename = os.path.join(StreamListDir, run)

            if os.path.exists(outputfilename) and os.path.exists(
                outputinfoname
            ):
                continue

            cshfile = scriptDir + "run" + runN + ".csh"
            # print cshfile
            FILE = open(cshfile, "w")
            FILE.writelines("#!/bin/tcsh \n")
            FILE.writelines("source ~lijj16/mroot.csh \n")
            FILE.writelines(f"./AnalysisLSAlphas {inputfilename} {outputfilename} {outputinfoname}")
            FILE.close()
            os.system("chmod u+x " + cshfile)
            # os.system("hep_sub -g dyw -wt mid " + cshfile)
            # Add the current job script to job list
            cshfiles = cshfiles + " " + cshfile
            job_count = job_count + 1

    # This part is to create a parent script to revoke all job scripts by procid
    ensure_dir(basedir + "parent/" + site)
    parent_script = basedir + "parent/" + site + "/parent_script.sh"
    
    FILE = open(parent_script, "w")
    FILE.writelines("#!/bin/bash \n")
    FILE.writelines("job_scripts=(" + cshfiles + ")\n")
    FILE.writelines("${job_scripts[$1]}")
    FILE.close()

    os.system("chmod u+x " + parent_script)

    os.system(
        "hep_sub -g dyw "
        + parent_script
        + " -argu %{ProcId} -n "
        + str(job_count)
    )
    """
    for longer job: -wt mid
    """
    print(
        "hep_sub -g dyw -wt mid "
        + parent_script
        + " -argu %{ProcId} -n "
        + str(job_count)
    )
