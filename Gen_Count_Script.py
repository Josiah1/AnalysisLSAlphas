#!/usr/bin/env python3
"""
   Submit Ostw jobs for all the list files in the specified directory.
   Need to manually edit "site" and "StreamListDir".

   Feb. 2012
   Zhe Wang
"""

import sys, os
import math

if __name__ == "__main__":

    site = sys.argv[1]
    job_limit = math.inf
    try:
        job_limit = int(sys.argv[2])
    except:
        pass

    data = "H"

    # AppDir = "/dybfs2/users/sunhaozhe/TWin/Retw_lijj/src_flasher_nucorr/"
    # AppDir = "/dybfs2/users/sunhaozhe/nHBackEndAnalysis/SpallationN/TWin/Retw/src_retwSN/"
    AppDir = "/dybfs2/users/lijj16/TWin/Retw_lijj/src_FD"

    ListOfStreamListDir = [
        "/dybfs2/rec/P17B/goodrunlist_p17b_v4/good_v4_sync.ihep/" + site + "/",
        "/dybfs2/rec/P19A/goodrunlist_p19a_v5/good_v5_sync.ihep/" + site + "/",
        "/dybfs2/rec/P20A/goodrunlist_p20a_v3/good_v3_sync.ihep/" + site + "/",
        "/dybfs2/rec/P21A/goodrunlist_p21a_v1/good_v1_sync.ihep/" + site + "/",
    ]

    basedir = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_IBD_Extended/"
    #basedir = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_IBD_0p7MeV/"

    print(AppDir)
    print(basedir)

    cshfiles = ""
    job_count = 0

    flag_exceed = False

    for StreamListDir in ListOfStreamListDir:
        ListOfList = os.listdir(StreamListDir)
        for runIdx, run in enumerate(sorted(ListOfList)):
            ext = run[-5:]
            if ext == ".list":

                if job_count + 1 > job_limit:
                    flag_exceed = True
                    break

                base = run[:-5]
                runN = run[3:8]
                scriptDir = basedir + "Script/" + site + "/"
                outputDir = basedir + "Output/" + site + "/"
                logDir = basedir + "Log/" + site + "/"
                infoDir = basedir + "RunInfo/" + site + "/"

                if os.path.exists(outputDir + runN + ".TWin.root") and os.path.exists(
                    infoDir + runN + "RunInfo.txt"
                ):
                    continue

                cshfile = scriptDir + "run" + runN + ".csh"
                # print cshfile
                FILE = open(cshfile, "w")
                FILE.writelines("#!/bin/tcsh \n")
                FILE.writelines("source ~lijj16/mroot.csh \n")
                FILE.writelines("cd " + AppDir + " \n")

                FILE.writelines(
                    "./Retw "
                    + StreamListDir
                    + base
                    + ".list"
                    + "  "
                    + outputDir
                    + runN
                    + ".TWin.root "
                    + infoDir
                    + runN
                    + "RunInfo.txt "
                    + data
                    + " > "
                    + logDir
                    + runN
                    + ".log\n"
                )

                FILE.close()
                os.system("chmod u+x " + cshfile)
                # os.system("hep_sub -g dyw -wt mid " + cshfile)
                # Add the current job script to job list
                cshfiles = cshfiles + " " + cshfile
                job_count = job_count + 1
        if flag_exceed:
            break

    # This part is to create a parent script to revoke all job scripts by procid
    parent_script = basedir + "parent/" + site + "/parent_script.sh"
    FILE = open(parent_script, "w")
    FILE.writelines("#!/bin/bash \n")
    FILE.writelines("job_scripts=(" + cshfiles + ")\n")
    FILE.writelines("${job_scripts[$1]}")
    FILE.close()

    os.system("chmod u+x " + parent_script)

    os.system(
        "hep_sub -g dyw -wt mid "
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
