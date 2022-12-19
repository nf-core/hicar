#!/usr/bin/env python
import struct, io, os, joblib, argparse, sys
import numpy as np
from sklearn.isotonic import IsotonicRegression
from scipy import sparse
from scipy import stats
from scipy.ndimage import gaussian_filter
from numba import njit
from collections import defaultdict
import cooler

# file copied data: https://github.com/tariks/peakachu/commit/adc736627bf43451aa1eee8ece061f3d57bc0c64
def tocsr(X):

    row, col, data = X.row, X.col, X.data
    M = sparse.csr_matrix((data, (row, col)), shape=X.shape, dtype=float)

    return M


def calculate_expected(M, maxdis, raw=False):

    n = M.shape[0]
    R, C = M.nonzero()
    valid_pixels = np.isfinite(M.data)
    # extract valid columns
    if raw:
        R, C, data = R[valid_pixels], C[valid_pixels], M.data[valid_pixels]
        M = sparse.csr_matrix((data, (R, C)), shape=M.shape, dtype=float)
        marg = np.array(M.sum(axis=0)).ravel()
        valid_cols = marg > 0
    else:
        R, C = set(R[valid_pixels]), set(C[valid_pixels])
        valid_cols = np.zeros(n, dtype=bool)
        for i in R:
            valid_cols[i] = True
        for i in C:
            valid_cols[i] = True

    # calculate the expected value for each genomic distance
    exp_arr = np.zeros(maxdis + 1)
    for i in range(maxdis + 1):
        if i == 0:
            valid = valid_cols
        else:
            valid = valid_cols[:-i] * valid_cols[i:]

        diag = M.diagonal(i)
        diag = diag[valid]
        if diag.size > 10:
            exp = diag.mean()
            exp_arr[i] = exp

    # make exp_arr stringently non-increasing
    IR = IsotonicRegression(increasing=False, out_of_bounds="clip")
    _d = np.where(exp_arr > 0)[0]
    IR.fit(_d, exp_arr[_d])
    exp_arr = IR.predict(list(range(maxdis + 1)))

    return exp_arr


@njit
def distance_normaize_core(sub, exp_bychrom, x, y, w):

    # calculate x and y indices
    x_arr = np.arange(x - w, x + w + 1).reshape((2 * w + 1, 1))
    y_arr = np.arange(y - w, y + w + 1)

    D = y_arr - x_arr
    D = np.abs(D)
    min_dis = D.min()
    max_dis = D.max()
    if max_dis >= exp_bychrom.size:
        return sub
    else:
        exp_sub = np.zeros(sub.shape)
        for d in range(min_dis, max_dis + 1):
            xi, yi = np.where(D == d)
            for i, j in zip(xi, yi):
                exp_sub[i, j] = exp_bychrom[d]

        normed = sub / exp_sub

        return normed


@njit
def image_normalize(arr_2d):

    arr_2d = (arr_2d - arr_2d.min()) / (arr_2d.max() - arr_2d.min())  # value range: [0,1]

    return arr_2d


@njit
def distance_normalize(arr_pool, exp_bychrom, xi, yi, w):

    clist = []
    fea = []
    for i in range(xi.size):
        x = xi[i]
        y = yi[i]
        window = arr_pool[i]

        bad_x, bad_y = np.where(np.isnan(window))
        for i_, j_ in zip(bad_x, bad_y):
            window[i_, j_] = 0

        if np.count_nonzero(window) < window.size * 0.1:
            continue

        ll_mean = window[:w, :w].mean()
        if ll_mean > 0:
            center = window[w, w]
            p2LL = center / ll_mean
            if p2LL > 0.1:
                window = distance_normaize_core(window, exp_bychrom, x, y, w)
                fea.append(window)
                clist.append((x, y))

    return fea, clist


class Chromosome:
    def __init__(self, M, model, raw_M=None, weights=None, lower=6, upper=300, cname="chrm", res=10000, width=5):

        lower = max(lower, width + 1)
        upper = min(upper, M.shape[0] - 2 * width)
        # calculate expected values
        if weights is None:
            self.exp_arr = calculate_expected(M, upper + 2 * width, raw=True)
            if M is raw_M:
                self.background = self.exp_arr
            else:
                self.background = calculate_expected(raw_M, upper + 2 * width, raw=True)
        else:
            self.exp_arr = calculate_expected(M, upper + 2 * width, raw=False)
            self.background = self.exp_arr

        self.raw_M = raw_M
        self.weights = weights

        # lower down the memory usage
        R, C = M.nonzero()
        validmask = np.isfinite(M.data) & (C - R > (-2 * width)) & (C - R < (upper + 2 * width))
        R, C, data = R[validmask], C[validmask], M.data[validmask]
        self.M = sparse.csr_matrix((data, (R, C)), shape=M.shape)
        self.get_candidate(lower, upper)
        self.chromname = cname
        self.r = res
        self.w = width
        self.model = model

    def get_candidate(self, lower, upper):

        x_arr = np.array([], dtype=int)
        y_arr = np.array([], dtype=int)
        p_arr = np.array([], dtype=float)
        idx = np.arange(self.raw_M.shape[0])
        for i in range(lower, upper + 1):
            diag = self.raw_M.diagonal(i)
            e = self.background[i]
            if (diag.size > 0) and (e > 0):
                xi = idx[:-i]
                yi = idx[i:]
                if self.weights is None:
                    exp = np.ones(diag.size, dtype=float) * e
                else:
                    b1 = self.weights[:-i]
                    b2 = self.weights[i:]
                    exp = np.ones(diag.size, dtype=float) * e / (b1 * b2)

                Poiss = stats.poisson(exp)
                pvalues = Poiss.sf(diag)
                mask = (diag > 0) & np.isfinite(pvalues)
                x_arr = np.r_[x_arr, xi[mask]]
                y_arr = np.r_[y_arr, yi[mask]]
                p_arr = np.r_[p_arr, pvalues[mask]]

        # qvalues = multipletests(p_arr, method = 'fdr_bh')[1]
        mask = p_arr < 0.01
        self.ridx, self.cidx = x_arr[mask], y_arr[mask]

    def getwindow(self, coords):

        w = self.w
        coords = np.r_[coords]
        xi, yi = coords[:, 0], coords[:, 1]
        mask = (xi - w >= 0) & (yi + w + 1 <= self.M.shape[0])
        xi, yi = xi[mask], yi[mask]
        seed = np.arange(-w, w + 1)
        delta = np.tile(seed, (seed.size, 1))
        xxx = xi.reshape((xi.size, 1, 1)) + delta.T
        yyy = yi.reshape((yi.size, 1, 1)) + delta
        v = np.array(self.M[xxx.ravel(), yyy.ravel()]).ravel()
        vvv = v.reshape((xi.size, seed.size, seed.size))
        windows, clist = distance_normalize(vvv, self.exp_arr, xi, yi, w)
        fea = []
        for arr in windows:
            tmp = gaussian_filter(arr, sigma=1, order=0)
            scaled_arr = image_normalize(tmp)
            fea.append(scaled_arr.ravel())

        fea = np.r_[fea]
        clist = np.r_[clist]

        return fea, clist

    def score(self, thre=0.5):

        print("scoring matrix {}".format(self.chromname))
        print("number of candidates {}".format(self.ridx.size))
        total_coords = [(r, c) for r, c in zip(self.ridx, self.cidx)]
        ri = np.r_[[]]
        ci = np.r_[[]]
        prob_pool = np.r_[[]]
        # to lower down the memory usage
        batch_size = 100000
        for t in range(0, len(total_coords), batch_size):
            coords = total_coords[t : t + batch_size]
            fea, clist = self.getwindow(coords)
            if fea.shape[0] > 1:
                p = self.model.predict_proba(fea)[:, 1]
                pfilter = p > thre
                ri = np.r_[ri, clist[:, 0][pfilter]]
                ci = np.r_[ci, clist[:, 1][pfilter]]
                prob_pool = np.r_[prob_pool, p[pfilter]]
        ri = ri.astype(int)
        ci = ci.astype(int)
        """
        # finely tune the probability score by combining the Fold-enrichment score
        prob_tuned = []
        for i in range(ri.size):
            fc = fc_pool[i]
            # change 0.1 and 0.9 if more weithts need to be given to FC
            fc_factor = (1 / (1 + np.exp(-fc)) - 0.5) / 0.5 * 0.1 + 0.9
            p = prob_pool[i] * fc_factor
            prob_tuned.append(p)
        prob_pool = np.r_[prob_tuned]
        """
        result = sparse.csr_matrix((prob_pool, (ri, ci)), shape=self.M.shape)
        data = np.array(self.M[ri, ci]).ravel()
        self.M = sparse.csr_matrix((data, (ri, ci)), shape=self.M.shape)

        return result, self.M

    def writeBed(self, outfil, prob_csr, raw_csr):

        with open(outfil, "a") as out:
            r, c = prob_csr.nonzero()
            for i in range(r.size):
                line = [
                    self.chromname,
                    r[i] * self.r,
                    (r[i] + 1) * self.r,
                    self.chromname,
                    c[i] * self.r,
                    (c[i] + 1) * self.r,
                    ".",
                    prob_csr[r[i], c[i]],
                    ".",
                    ".",
                    raw_csr[r[i], c[i]],
                ]
                out.write("\t".join(list(map(str, line))) + "\n")


def getargs():
    ## Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(
        description="""Unveil Hi-C Anchors and Peaks.""", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    for i in subs[1:]:
        i.add_argument("-r", "--resolution", help="Resolution in bp (default 10000)", type=int, default=10000)
    for i in subs[:-1]:
        i.add_argument("-p", "--path", help="Path to a .cool URI string or a .hic file.")

    for i in subs[1:-1]:
        i.add_argument("--balance", action="store_true", help="""Whether or not using the ICE/KR-balanced matrix.""")

    subchrom.add_argument(
        "-C",
        "--chrom",
        help="""Chromosome label. Only contact data within the
                          specified chromosome will be considered.""",
    )
    subgen.add_argument(
        "-C",
        "--chroms",
        nargs="*",
        default=["#", "X"],
        help="List of chromosome labels. Only contact data within the specified "
        'chromosomes will be included. Specially, "#" stands for chromosomes '
        'with numerical labels. "--chroms" with zero argument will include '
        'all chromosome data. (default "#" X)',
    )

    for i in subs[2:-1]:
        i.add_argument("-m", "--model", type=str, help="""Path to pickled model file.""")
        i.add_argument(
            "-l", "--lower", type=int, default=6, help="""Lower bound of distance between loci in bins (default 6)."""
        )
        i.add_argument(
            "-u",
            "--upper",
            type=int,
            default=300,
            help="""Upper bound of distance between loci in bins (default 300).""",
        )
        i.add_argument(
            "--minimum-prob",
            type=float,
            default=0.5,
            help="""Only output pixels with probability score greater than this value (default 0.5)""",
        )
        i.add_argument("-O", "--output", help="Output file name.")

    subdepth.add_argument(
        "--min-dis",
        default=0,
        type=int,
        help="""Only count reads with genomic distance (in base pairs) greater than this value. (default 0)""",
    )

    subtrain.add_argument("-b", "--bedpe", help="""Path to the bedpe file containing positive training set.""")
    subtrain.add_argument(
        "-w",
        "--width",
        type=int,
        default=5,
        help="""Number of bins added to center of window.
                        default width=5 corresponds to 11x11 windows""",
    )
    subtrain.add_argument(
        "--nproc",
        type=int,
        default=4,
        help="""Number of worker processes that
                          will be allocated for training. (default 4)""",
    )
    subtrain.add_argument("-O", "--output", help="Folder path to store trained models.")

    subpool.add_argument(
        "-i", "--infile", help="""Path to the bedpe file outputted from score_chromosome or score_genome"""
    )
    subpool.add_argument("-o", "--outfile", help="Output file name.")
    subpool.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.9,
        help="Probability threshold applied before peak calling (default 0.9)",
    )

    ## Parse the command-line arguments
    commands = sys.argv[1:]
    if (not commands) or (
        (commands[0] in ["train", "score_chromosome", "score_genome", "depth", "pool"]) and len(commands) == 1
    ):
        commands.append("-h")
    args = parser.parse_args(commands)

    return args, commands


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--resolution", help="Resolution in bp (default 10000)", type=int, default=10000)
    parser.add_argument("-p", "--path", help="Path to a .cool URI string or file.")
    parser.add_argument("--balance", action="store_true", help="""Whether or not using the ICE/KR-balanced matrix.""")
    parser.add_argument(
        "-C",
        "--chroms",
        nargs="*",
        default=["#", "X"],
        help="List of chromosome labels. Only contact data within the specified "
        'chromosomes will be included. Specially, "#" stands for chromosomes '
        'with numerical labels. "--chroms" with zero argument will include '
        'all chromosome data. (default "#" X)',
    )
    parser.add_argument("-m", "--model", type=str, help="""Path to pickled model file.""")
    parser.add_argument(
        "-l", "--lower", type=int, default=6, help="""Lower bound of distance between loci in bins (default 6)."""
    )
    parser.add_argument(
        "-u", "--upper", type=int, default=300, help="""Upper bound of distance between loci in bins (default 300)."""
    )
    parser.add_argument(
        "--minimum-prob",
        type=float,
        default=0.5,
        help="""Only output pixels with probability score greater than this value (default 0.5)""",
    )
    parser.add_argument("-O", "--output", help="Output file name.")

    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    np.seterr(divide="ignore", invalid="ignore")

    if os.path.exists(args.output):
        os.remove(args.output)

    model = joblib.load(args.model)

    # deduce the width parameter used during the training
    width = int((np.sqrt(model.feature_importances_.size) - 1) / 2)

    # not support .hic
    Lib = cooler.Cooler(args.path)
    chromosomes = Lib.chromnames[:]

    queue = []
    for key in chromosomes:

        chromlabel = key.lstrip("chr")
        if (not args.chroms) or (chromlabel.isdigit() and "#" in args.chroms) or (chromlabel in args.chroms):
            queue.append(key)

    for key in queue:

        if key.startswith("chr"):
            cname = key
        else:
            cname = "chr" + key

        if args.balance:
            M = tocsr(Lib.matrix(balance=args.balance, sparse=True).fetch(key))
            raw_M = tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
            weights = Lib.bins().fetch(key)["weight"].values
            X = Chromosome(
                M,
                model=model,
                raw_M=raw_M,
                weights=weights,
                cname=cname,
                lower=args.lower,
                upper=args.upper,
                res=args.resolution,
                width=width,
            )
        else:
            M = tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
            X = Chromosome(
                M,
                model=model,
                raw_M=M,
                weights=None,
                cname=cname,
                lower=args.lower,
                upper=args.upper,
                res=args.resolution,
                width=width,
            )

        result, R = X.score(thre=args.minimum_prob)
        X.writeBed(args.output, result, R)
        return 0


if __name__ == "__main__":
    sys.exit(main())
