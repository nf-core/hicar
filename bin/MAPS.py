#!/usr/bin/env python

#########################################
# Author: [Ivan Juric](https://github.com/ijuric)
# File: MAPS.py
# Source: https://github.com/ijuric/MAPS/blob/master/bin/MAPS/MAPS.py
# Source+commit: https://raw.githubusercontent.com/ijuric/MAPS/860f864401b31a085f2a091ba593cf579442ac35/bin/MAPS/MAPS.py
# Data: 11/08/2021, commit: 860f864
# Copyright (c)  Ivan Juric ( jurici@ccf.org ) and Ming Hu ( hum@ccf.org )
# This source code is licensed under the GPL-3.0 license
# modified by Jianhong:
# ## 1. reset the python environment.
# ## 2. Handle the error when MACS2 input is empty
# ## 3. Handle the error when the metadata file is empty
# ## 4. Automatic detect the chromosome names, function: get_chrom_from_MACS2
# ## 5. Fix the chromosome style issues
# ## 6. Handle the NA values
# ## 7. Fix the wrong data type when reading chromosome name as number
# ## 8. add weight for counts in case the input file is a count ready file
# ## 9. fix the indent space
# ##10. add count for inter-chromosomes, and fix the filenames
#########################################

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
import sys
import itertools
import re


def get_segment_range(bin_start, bin_end):
    return range(int(bin_start), int(bin_end) + 1)


def get_1D_peak(segment_range, MACS2_peak_ranges_list):
    for segment in segment_range:
        if segment in MACS2_peak_ranges_list:
            return True
    return False


def validate_input_data(input_data, long_bedpe_postfix, short_bed_postfix):
    params = {}
    if "OUT_DIR" in input_data:
        params["OUT_DIR"] = input_data["OUT_DIR"][1]
    if "BINNING_RANGE" in input_data:
        params["BINNING_RANGE"] = input_data["BINNING_RANGE"].astype("int")[1]
    if "BIN_SIZE" in input_data:
        params["BIN_SIZE"] = input_data["BIN_SIZE"].astype("int")[1]
    if "DATASET_NAME" in input_data:
        params["DATASET_NAME"] = input_data["DATASET_NAME"][1]
    if "MACS2_PATH" in input_data:
        params["MACS2_PATH"] = input_data["MACS2_PATH"][1]
    if "GF_PATH" in input_data:
        params["GF_PATH"] = input_data["GF_PATH"][1]
    if "LONG_PATH" in input_data:
        params["LONG_PATH"] = input_data["LONG_PATH"][1]
    if "LONG_FORMAT" in input_data:
        params["LONG_FORMAT"] = input_data["LONG_FORMAT"][1]
        if "long.intra.bedpe" in params["LONG_FORMAT"]:
            params["LONG_FORMAT"] = params["LONG_FORMAT"].replace("long.intra.bedpe", long_bedpe_postfix)
    if "SHORT_PATH" in input_data:
        params["SHORT_PATH"] = input_data["SHORT_PATH"][1]
    if "SHORT_FORMAT" in input_data:
        params["SHORT_FORMAT"] = input_data["SHORT_FORMAT"][1]
        if "shrt.vip.bed" in params["SHORT_FORMAT"]:
            params["SHORT_FORMAT"] = params["SHORT_FORMAT"].replace("shrt.vip.bed", short_bed_postfix)
    # if 'N_CHROMS' in input_data:
    #    params['N_CHROMS'] = input_data['N_CHROMS'].astype('int')[1]
    if "SEX_CHROMS" in input_data:
        params["SEX_CHROMS"] = input_data["SEX_CHROMS"][1]
    else:
        params["SEX_CHROMS"] = ""
    return params


def pd_read_tab(*coln, **args):
    try:
        out = pd.read_csv(**args)
    except Exception:
        out = pd.DataFrame(columns=coln)
    return out


def load_MACS2(MACS2_PATH):
    coln = ["chr", "start", "end"]
    MACS2_full = pd_read_tab(
        coln,
        filepath_or_buffer=MACS2_PATH,
        sep="\t",
        skip_blank_lines=True,
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        low_memory=False,
    )
    MACS2_full.columns = coln
    MACS2_full = MACS2_full.astype({"chr": str})
    return MACS2_full


def load_metadata(GF_PATH, BIN_SIZE):
    coln = ["chr", "start", "end", "effective_length", "gc", "mappability", "bin1_mid", "bin2_mid", "bin"]
    metadata_full = pd_read_tab(coln, filepath_or_buffer=GF_PATH, sep="\t", header=None, low_memory=False)
    metadata_full.columns = ["chr", "start", "end", "effective_length", "gc", "mappability"]
    metadata_full = metadata_full.astype({"chr": str})
    metadata_full["bin1_mid"] = metadata_full["start"] // BIN_SIZE
    metadata_full["bin2_mid"] = metadata_full["bin1_mid"]
    metadata_full["bin"] = metadata_full["bin1_mid"]
    return metadata_full


def parse_fname(chrom, type, params):
    if type == "short":
        fname = params["SHORT_PATH"] + params["SHORT_FORMAT"]
    else:
        fname = params["LONG_PATH"] + params["LONG_FORMAT"]
    for p in params:
        pn = "[" + p + "]"
        if pn in fname:
            fname = fname.replace(pn, params[p])
    if "[CHROMOSOME]" in fname:
        fname = fname.replace("[CHROMOSOME]", chrom)
    else:
        if type == "short":
            print("File format needs to contain [CHROMOSOME] tag:", params["SHORT_FORMAT"])
        else:
            print("File format needs to contain [CHROMOSOME] tag:", params["LONG_FORMAT"])
        exit()
    return fname


def get_chrom_from_MACS2(MACS2_full):
    chr = []
    if MACS2_full.shape[0]:
        chr = MACS2_full["chr"].unique().tolist()
        ## remove XY, MT, ...
        p = re.compile("^(chr)?((\d{1,2})|(IX|IV|V?I{0,3}))$")
        chr = [s for s in chr if p.match(s)]
    return chr


def get_peaks_range(MACS2_full, CHR, params):
    MACS2 = MACS2_full[MACS2_full["chr"] == CHR].copy()
    if not MACS2.shape[0]:
        CHR = CHR.replace("chr", "")
        MACS2 = MACS2_full[MACS2_full["chr"] == CHR].copy()
        if not MACS2.shape[0]:
            return True, False
    MACS2["start_bin"] = np.floor(MACS2["start"] / params["BIN_SIZE"]).fillna(0)
    MACS2["end_bin"] = np.ceil(MACS2["end"] / params["BIN_SIZE"]).fillna(0)
    # perform this hack becasue apply returns wrong data type in some rare case
    specialCase = False
    if MACS2.iloc[0]["end_bin"] - MACS2.iloc[0]["start_bin"] == MACS2.shape[1] - 1:
        MACS2.iloc[0, MACS2.columns.get_loc("start_bin")] = MACS2.iloc[0]["start_bin"] - 1
        specialCase = True
    MACS2_peak_ranges = MACS2.apply(
        lambda row: range(int(row["start_bin"]), int(row["end_bin"])), axis=1
    ).values.tolist()
    MACS2_peak_ranges_list = set(itertools.chain.from_iterable(MACS2_peak_ranges))
    if specialCase:
        MACS2_peak_ranges_list.remove(MACS2.iloc[0]["start_bin"])
    return False, MACS2_peak_ranges_list


def init(p):
    ## checking that all files are available
    print("loading parameters file")
    input_data = pd_read_tab(
        ["chr", "start", "end"],
        filepath_or_buffer=p.run_file,
        sep="=",
        skip_blank_lines=True,
        comment="#",
        index_col=0,
        header=None,
    )
    input_data = input_data.transpose()
    params = validate_input_data(input_data, p.long_bedpe_postfix, p.short_bed_postfix)
    print("loading MACS2 peaks")
    MACS2_full = load_MACS2(params["MACS2_PATH"])
    ## setting up chromosomes, TODO: extract the chromosome from MACS2_full["chr"], make sure use the clean chromosomes
    # chroms = ['chr' + str(i) for i in range(1,params['N_CHROMS']+1,1)]
    chroms = get_chrom_from_MACS2(MACS2_full)
    print(chroms)
    if params["SEX_CHROMS"] == "X":
        chroms.append("chrX")
    elif params["SEX_CHROMS"] == "Y":
        chroms.append("chrY")
    elif params["SEX_CHROMS"] == "XY":
        chroms.append("chrX")
        chroms.append("chrY")
    print(chroms)
    params["BIN_RANGE"] = float(params["BINNING_RANGE"]) / float(params["BIN_SIZE"])
    print("loading metadata file")
    metadata_full = load_metadata(params["GF_PATH"], params["BIN_SIZE"])
    qc_str = ""  ## content of qc.maps file
    for CHR1 in chroms:
        peak_skip, MACS2_peak_ranges_list_1 = get_peaks_range(MACS2_full, CHR1, params)
        ps_short1 = pd_read_tab(
            ["chr", "start", "end"], filepath_or_buffer=parse_fname(CHR1, "short", params), header=None, sep="\t"
        )
        if peak_skip:
            continue
        for CHR2 in chroms:
            print("doing chromosome ", CHR1, " and ", CHR2, "\n")
            # handling MACS2 peaks
            print("-- handling MACS2 peaks")
            peak_skip2, MACS2_peak_ranges_list_2 = get_peaks_range(MACS2_full, CHR2, params)
            if peak_skip2:
                continue
            print("-- handling short.bed\n")
            if CHR1 != CHR2:
                ps_short2 = pd_read_tab(
                    ["chr", "start", "end"],
                    filepath_or_buffer=parse_fname(CHR2, "short", params),
                    header=None,
                    sep="\t",
                )
                ps_short = pd.concat([ps_short1, ps_short2], ignore_index=True)
            else:
                ps_short = ps_short1
            if ps_short.shape[0]:
                new_cols = ["chr", "start", "end", "name"]
                ps_short.rename(columns=dict(zip(ps_short.columns[0:], new_cols)), inplace=True)
                ps_short = ps_short.astype({"chr": str})
                ps_short["bin"] = ps_short[["start", "end"]].mean(axis=1) // params["BIN_SIZE"]
                ps_short["short_count"] = 1
                count_data_short = ps_short[["chr", "bin", "short_count"]].groupby(["chr", "bin"]).count()
                count_data_short.reset_index(inplace=True)
                print("-- handling long.bedpe\n")
                ##### getting overlap
                ## load long.bed file
                long_cols = ["chr1", "start1", "end1", "chr2", "start2", "end2", "count"]
                ps_long = pd_read_tab(
                    long_cols,
                    filepath_or_buffer=parse_fname(CHR1 + "_" + CHR2, "long", params),
                    header=None,
                    sep="\t",
                    low_memory=False,
                )
                ps_long.rename(columns=dict(zip(ps_long.columns[0:], long_cols)), inplace=True)
                if ps_long.shape[0]:
                    ps_long = ps_long.astype({"chr1": str, "chr2": str})
                    ## filter only reads at the same chromosome and proper orientation
                    ps_long = ps_long[(ps_long["chr1"] == CHR1) & (ps_long["chr2"] == CHR2)]
                if ps_long.shape[0]:
                    ps_long["read1_bin_mid"] = ((ps_long["start1"] + ps_long["end1"]) / 2.0) // params["BIN_SIZE"]
                    ps_long["read2_bin_mid"] = ((ps_long["start2"] + ps_long["end2"]) / 2.0) // params["BIN_SIZE"]
                    ps_long["bin1_mid"] = ps_long.loc[:, ["read1_bin_mid", "read2_bin_mid"]].min(axis=1)
                    ps_long["bin2_mid"] = ps_long.loc[:, ["read1_bin_mid", "read2_bin_mid"]].max(axis=1)
                    # ps_long['count'] = 1
                    # count_data = ps_long[['bin1_mid', 'bin2_mid','count']].groupby(['bin1_mid','bin2_mid']).count()
                    count_data = ps_long[["bin1_mid", "bin2_mid", "count"]]
                    count_data.reset_index(inplace=True)
                    count_data_and = count_data[
                        (count_data["bin1_mid"].isin(MACS2_peak_ranges_list_1))
                        & (count_data["bin2_mid"].isin(MACS2_peak_ranges_list_2))
                    ].copy()
                    if CHR1 == CHR2:
                        count_data_and = count_data_and[
                            (np.abs(count_data_and["bin1_mid"] - count_data_and["bin2_mid"]) <= params["BIN_RANGE"])
                            & (np.abs(count_data_and["bin1_mid"] - count_data_and["bin2_mid"]) >= 1)
                        ]
                    count_data_and["1D_peak_bin1"] = 1
                    count_data_and["1D_peak_bin2"] = 1
                    count_data_xor = count_data[
                        (count_data.bin1_mid.isin(MACS2_peak_ranges_list_1))
                        ^ (count_data.bin2_mid.isin(MACS2_peak_ranges_list_2))
                    ]
                    if CHR1 == CHR2:
                        count_data_xor = count_data_xor[
                            (np.abs(count_data_xor["bin1_mid"] - count_data_xor["bin2_mid"]) <= params["BIN_RANGE"])
                            & (np.abs(count_data_xor["bin1_mid"] - count_data_xor["bin2_mid"]) >= 1)
                        ]
                    count_data_xor_bin1 = count_data_xor[
                        (count_data_xor.bin1_mid.isin(MACS2_peak_ranges_list_1))
                    ].copy()
                    count_data_xor_bin1["1D_peak_bin1"] = 1
                    count_data_xor_bin1["1D_peak_bin2"] = 0
                    count_data_xor_bin2 = count_data_xor[
                        (count_data_xor.bin2_mid.isin(MACS2_peak_ranges_list_2))
                    ].copy()
                    count_data_xor_bin2["1D_peak_bin1"] = 0
                    count_data_xor_bin2["1D_peak_bin2"] = 1
                    count_data_xor = pd.concat([count_data_xor_bin1, count_data_xor_bin2], ignore_index=True)
                    print("-- calculating values for maps.qc file\n")
                    AND_sum = count_data_and["count"].sum()
                    XOR_sum = count_data_xor["count"].sum()
                    NOT_sum = count_data["count"].sum() - AND_sum - XOR_sum
                    qc_str = (
                        qc_str
                        + "AND_set\t"
                        + str(AND_sum)
                        + "\tnumber of pairs in AND set at chromsome "
                        + CHR1
                        + " and "
                        + CHR2
                        + "\n"
                        + "XOR_set\t"
                        + str(XOR_sum)
                        + "\tnumber of pairs in XOR set at chromsome "
                        + CHR1
                        + " and "
                        + CHR2
                        + "\n"
                        + "NOT_set\t"
                        + str(NOT_sum)
                        + "\tnumber of pairs in NOT set at chromsome "
                        + CHR1
                        + " and "
                        + CHR2
                        + "\n"
                    )
                    print("-- handling metadata\n")
                    metadata = metadata_full[metadata_full["chr"].isin([CHR1, CHR2])].copy()
                    metadata = pd.merge(metadata, count_data_short, on=["bin", "chr"], how="outer")
                    metadata["short_count"] = metadata["short_count"].fillna(0)
                    print("-- attaching genome features atributes to AND set")
                    reg_and = pd.merge(
                        count_data_and,
                        metadata[["bin1_mid", "effective_length", "gc", "mappability", "short_count"]],
                        on="bin1_mid",
                    )
                    reg_and.rename(
                        columns={
                            "effective_length": "effective_length1",
                            "gc": "gc1",
                            "mappability": "mappability1",
                            "short_count": "short_count1",
                        },
                        inplace=True,
                    )
                    reg_and = pd.merge(
                        reg_and,
                        metadata[["bin2_mid", "effective_length", "gc", "mappability", "short_count"]],
                        on="bin2_mid",
                    )
                    reg_and.rename(
                        columns={
                            "effective_length": "effective_length2",
                            "gc": "gc2",
                            "mappability": "mappability2",
                            "short_count": "short_count2",
                        },
                        inplace=True,
                    )
                    reg_and = reg_and[(reg_and["effective_length1"] > 0) & (reg_and["effective_length2"] > 0)]
                    if CHR1 == CHR2:
                        reg_and["dist"] = pd.to_numeric(np.abs(reg_and["bin1_mid"] - reg_and["bin2_mid"]))
                    else:
                        reg_and["dist"] = 9223372036854775807
                    reg_and["logl"] = np.log(
                        (reg_and["effective_length1"] + 1.0)
                        * (reg_and["effective_length2"] + 1.0)
                        / (params["BIN_SIZE"] * params["BIN_SIZE"])
                    )
                    reg_and["loggc"] = np.log(reg_and["gc1"] * reg_and["gc2"])
                    reg_and["logm"] = np.log(reg_and["mappability1"] * reg_and["mappability2"])
                    reg_and["logdist"] = np.log((1.0 + reg_and["dist"]) / params["BIN_RANGE"])
                    max_short_and = (reg_and["short_count1"].max() + 1.0) * (reg_and["short_count2"].max() + 1.0)
                    reg_and["logShortCount"] = np.log(
                        (reg_and["short_count1"] + 1.0) * (reg_and["short_count2"] + 1.0) / max_short_and
                    )
                    reg_and["bin1_mid"] = reg_and["bin1_mid"] * params["BIN_SIZE"]
                    reg_and["bin2_mid"] = reg_and["bin2_mid"] * params["BIN_SIZE"]
                    print("-- attaching genome features atributes to XOR set")
                    reg_xor = pd.merge(
                        count_data_xor,
                        metadata[["bin1_mid", "effective_length", "gc", "mappability", "short_count"]],
                        on="bin1_mid",
                    )
                    reg_xor.rename(
                        columns={
                            "effective_length": "effective_length1",
                            "gc": "gc1",
                            "mappability": "mappability1",
                            "short_count": "short_count1",
                        },
                        inplace=True,
                    )
                    reg_xor = pd.merge(
                        reg_xor,
                        metadata[["bin2_mid", "effective_length", "gc", "mappability", "short_count"]],
                        on="bin2_mid",
                    )
                    reg_xor.rename(
                        columns={
                            "effective_length": "effective_length2",
                            "gc": "gc2",
                            "mappability": "mappability2",
                            "short_count": "short_count2",
                        },
                        inplace=True,
                    )
                    reg_xor = reg_xor[(reg_xor["effective_length1"] > 0) & (reg_xor["effective_length2"] > 0)]
                    reg_xor["dist"] = pd.to_numeric(np.abs(reg_xor["bin1_mid"] - reg_xor["bin2_mid"]))
                    reg_xor["logl"] = np.log(
                        (reg_xor["effective_length1"] + 1.0)
                        * (reg_xor["effective_length2"] + 1.0)
                        / (params["BIN_SIZE"] * params["BIN_SIZE"])
                    )
                    reg_xor["loggc"] = np.log(reg_xor["gc1"] * reg_xor["gc2"])
                    reg_xor["logm"] = np.log(reg_xor["mappability1"] * reg_xor["mappability2"])
                    reg_xor["logdist"] = np.log((1.0 + reg_xor["dist"]) / params["BIN_RANGE"])
                    max_short_xor = (reg_xor["short_count1"].max() + 1.0) * (reg_xor["short_count2"].max() + 1.0)
                    reg_xor["logShortCount"] = np.log(
                        (reg_xor["short_count1"] + 1.0) * (reg_xor["short_count2"] + 1.0) / max_short_xor
                    )
                    reg_xor["bin1_mid"] = reg_xor["bin1_mid"] * params["BIN_SIZE"]
                    reg_xor["bin2_mid"] = reg_xor["bin2_mid"] * params["BIN_SIZE"]
                    print("--saving output\n")
                    fout_name = (
                        params["OUT_DIR"]
                        + "reg_raw."
                        + str(CHR1)
                        + "_"
                        + str(CHR2)
                        + "."
                        + params["DATASET_NAME"]
                        + "."
                        + str(int(params["BIN_SIZE"] / 1000))
                        + "k.and"
                    )
                    reg_and.to_csv(fout_name, sep="\t")
                    fout_name = (
                        params["OUT_DIR"]
                        + "reg_raw."
                        + str(CHR1)
                        + "_"
                        + str(CHR2)
                        + "."
                        + params["DATASET_NAME"]
                        + "."
                        + str(int(params["BIN_SIZE"] / 1000))
                        + "k.xor"
                    )
                    reg_xor.to_csv(fout_name, sep="\t")
            else:
                print(
                    "no bin pairs in long or short bedpe files for chromosome ",
                    CHR1,
                    " and ",
                    CHR2,
                    ". Doing next chromosome",
                )
    print("-- saving .qc.maps file\n")
    qc_fname = params["OUT_DIR"] + params["DATASET_NAME"] + ".maps.qc"
    qc_file = open(qc_fname, "w")
    qc_file.write(qc_str)
    qc_file.close()


def main():
    parser = argparse.ArgumentParser()
    parser.prog = "PROG"
    parser.description = "MAPS"
    parser.epilog = "This is where the command-line utility's epilog goes."
    parser.add_argument("run_file", help="file containing run parameters")
    parser.add_argument("long_bedpe_postfix", help="file extension for long bedpe")
    parser.add_argument("short_bed_postfix", help="file extension for short bed")
    p = parser.parse_args(sys.argv[1:])
    init(p)


if __name__ == "__main__":
    main()
