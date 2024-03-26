import argparse as _argparse
import pymadx as _pymadx

def main_PrintRMatrix():
    parser = _argparse.ArgumentParser()
    parser.add_argument("tfsfile", type=str, help="TFS file name")
    args = parser.parse_args()
    d = _pymadx.Data.Tfs(args.tfsfile)
    _pymadx.Plot.PrintRMatrixTable(d)


def main_PrintRMatrixPdf():
    parser = _argparse.ArgumentParser()
    parser.add_argument("tfsfile", type=str, help="TFS file name")
    parser.add_argument("outputname", type=str, help="Output name for pdf")
    args = parser.parse_args()

    outname = str(args.outputname)
    if not outname.endswith(".pdf"):
        outname += ".pdf"

    tfs = _pymadx.Data.Tfs(args.tfsfile)
    _pymadx.Plot.RMatrixTableToPdf(tfs, outname)
