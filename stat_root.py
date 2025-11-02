import os
import sys
import ROOT

data1_DIR = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_IBD_0p7MeV/RunInfo/"
data2_DIR = "/dybfs2/users/lijj16/DATA/Retw/IBD/H_U235/RunInfo/"


def verify_timestamp(infodir, site, runNo):
    timestamp_cutingedge = 1388505600

    infofile = infodir + f"{site}/{runNo}RunInfo.txt"
    with open(infofile, "r") as f:
        timestamp = int(f.readline().split()[2])
        if timestamp > timestamp_cutingedge:
            return False
        else:
            return True


def read_live_time(infodir):
    live_time = {"EH1": [0, 0, 0, 0], "EH2": [0, 0, 0, 0], "EH3": [0, 0, 0, 0]}
    for site in sorted(os.listdir(infodir)):
        infodir_site = os.path.join(infodir, site)
        for runfile in sorted(os.listdir(infodir_site)):
            infofile_path = os.path.join(infodir_site, runfile)
            f = open(infofile_path, "r")
            line_content = f.readline().split()
            fulltime = float(line_content[6])
            for idx in range(4):
                livetime1 = float(line_content[7 + idx])
                live_time[site][idx] += fulltime - livetime1
    return live_time


def inte_nevts_by_ZR2(site, iad, isotope):
    rootfile_name = None
    if isotope == "Bi214" or isotope == "Bi212":
        rootfile_name = f"counting_{site}.root"
    elif isotope == "Rn219" or isotope == "Po210":
        rootfile_name = f"countingLowE_{site}.root"
    else:
        sys.stderr.write("Error: Unknown isotope!\n")
        return -1

    rootfile_obj = ROOT.TFile(rootfile_name, "READ")

    hist_obj = None
    if not isotope == "Po210":
        hist_obj = rootfile_obj.Get(f"h2dZR2Sig_{isotope}_{iad}")
    else:
        hist_obj = rootfile_obj.Get(f"h2dZR2_{isotope}_{iad}")

    criteria_set = {
        "GdLS": {"R2": (0, 2.3), "Z": (-1.5, 1.5)},
        "Total": {"R2": (0, 5.0), "Z": (-3.0, 3.0)},
    }

    inte_num_GdLS = hist_obj.Integral(
        hist_obj.GetXaxis().FindBin(criteria_set["GdLS"]["R2"][0]),
        hist_obj.GetXaxis().FindBin(criteria_set["GdLS"]["R2"][1]),
        hist_obj.GetYaxis().FindBin(criteria_set["GdLS"]["Z"][0]),
        hist_obj.GetYaxis().FindBin(criteria_set["GdLS"]["Z"][1]),
    )
    inte_num_Total = hist_obj.Integral(
        hist_obj.GetXaxis().FindBin(criteria_set["Total"]["R2"][0]),
        hist_obj.GetXaxis().FindBin(criteria_set["Total"]["R2"][1]),
        hist_obj.GetYaxis().FindBin(criteria_set["Total"]["Z"][0]),
        hist_obj.GetYaxis().FindBin(criteria_set["Total"]["Z"][1]),
    )
    inte_num_LS = inte_num_Total - inte_num_GdLS

    rootfile_obj.Close()

    return inte_num_GdLS, inte_num_LS, inte_num_Total


if __name__ == "__main__":
    live_times = read_live_time(data1_DIR)
    print("Live times:", live_times)

    sites = ["EH1", "EH2", "EH3"]
    isotopes = ["Bi214", "Bi212", "Rn219", "Po210"]

    for site in sites:
        for iad in [1, 2, 3, 4]:
            for isotope in isotopes:
                counts = inte_nevts_by_ZR2(site, iad, isotope)
                rates = [count / live_times[site][iad - 1] for count in counts]
                if counts[2] > 0:
                    print(
                        f"Rates for {isotope} at {site}, IAD {iad}: {rates[0]:.2e}(GdLS) {rates[1]:.2e}(LS) [Hz]",
                    )
