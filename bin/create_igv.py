#!/usr/bin/env python

#########################################
# Author: Jianhong Ou
#########################################

from __future__ import print_function
import argparse
import os
import sys
import io
import errno
import glob

TrackType = {'bw':'bigwig', 'bb':'bigBed', 'bigbed':'bigBed', 'gff':'gff3',
        'narrowpeak':'narrowPeak', 'broadpeak':'broadPeak'}
types = {'bed':'annotation', 'gff3':'annotation', 'gtf':'annotation',
        'narrowPeak':'annotation', 'broadPeak':'annotation',
        'bigBed':'annotation', 'maf':'mut', 'bigwig':'wig',
        'bedpe':'interact', 'gwas':'gwas', 'seg':'seg',
        'bam':'alignment', 'cram':'alignment'}
def getFormat(f):
    extension = os.path.splitext(f)[1].replace(".", "").lower()
    if extension in TrackType.keys():
        extension = TrackType[extension]
    return extension

def create_igv(ListFile, Genome, PublishDir):
    prefix = ""
    if PublishDir:
        depth = PublishDir.count('/')+1
        prefix = '/'.join(['..']*depth) + '/'
    fileList = []
    fin = io.open(ListFile, mode="r")
    while True:
        line = fin.readline()
        if line:
            ifile = [x.strip() for x in line.strip().split('\t')]
            fmt = getFormat(ifile[1])
            index = ""
            if len(ifile)>2:
                index = '"indexURL":"' + ifile[2] + '",'
            fileList.append(f"""
            {{
                \"name\":\"{ifile[0]}\",
                \"url\":\"{prefix}{ifile[1]}\",
                \"type\":\"{types[fmt]}\",
                {index}
                \"format\":\"{fmt}\",
                \"height\":50
            }}
            """)
        else:
            break
            fin.close()

    html = f"""var options =
    {{
        \"genome\": \"{Genome}\",
        \"tracks\":[
    """
    html = html + ','.join(fileList) + """
    ]
    };
    """
    return html


def wrap_html(contents):
    header = """<!DOCTYPE html><html>
    <head>
        <script src="https://cdn.jsdelivr.net/npm/igv@2.10.1/dist/igv.min.js"></script>
    </head>
    <body>
    <div id='igv-div'></div>
    <script>
    var igvDiv = document.getElementById("igv-div");
    """
    footer = """
    igv.createBrowser(igvDiv, options)
        .then(function (browser) {
            console.log("Created IGV browser");
        })
    </script>
    </body>
    </html>
    """
    return header + contents + footer


def parse_args(args=None):
    Description = 'Create IGV file from a list of files.'
    Epilog = """Example usage: python create_igv.py <LIST_FILE> <GENOME> <PUBLISH_DIR>"""
    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    ## REQUIRED PARAMETERS
    argParser.add_argument('LIST_FILE', help="Tab-delimited file containing two columns i.e. samplename\tsignalfile. Header isnt required.")
    argParser.add_argument('GENOME', help="Full path to genome fasta file or shorthand for genome available in UCSC e.g. hg19.")
    argParser.add_argument('PUBLISH_DIR', help="Publish dir for igv. It will be used to calculate the relative folder from publish dir to outdir.")
    return argParser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    tracks = create_igv(
                    ListFile=args.LIST_FILE,
                    Genome=args.GENOME,
                    PublishDir=args.PUBLISH_DIR)
    html = wrap_html(tracks)
    out = open("index.html", "a")
    out.write(html)
    out.close()


if __name__ == "__main__":
    sys.exit(main())
