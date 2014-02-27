import sys
import os
from collections import defaultdict
from toolshed import reader, nopen

def faseq(fa, chrom, start, end, cache=[None]):
    """
    since bam is ordered by chrom, we read each chrom into memory in cache
    """
    if cache[0] is None or cache[0][0] != chrom:
        cache[0] = (chrom, "".join(x.strip() for i, x in
            enumerate(nopen("|samtools faidx %s %s" % (fa, chrom))) if i >
            0).upper())
    chrom, seq = cache[0]
    #return seq[start - 1: end]
    return seq[start - 2: end + 1]

def check_bias(bam, reference, n_check=2e6, min_map_q=20):

    fpath = "%s.bias.png" % os.path.splitext(bam)[0]
    opath = "%s.bias.txt" % os.path.splitext(bam)[0]

    # keys of meth_counts are the base-position in the read
    meth_counts = [ # first read, 2nd read
            defaultdict(lambda : {'C': 0, 'T': 0, 'G': 0, 'A': 0}),
            defaultdict(lambda : {'C': 0, 'T': 0, 'G': 0, 'A': 0}) ]

    n_to_check = sys.maxint if n_check is None else n_check
    fh = open(opath, "w")
    for toks in reader("|samtools view -F4 {bam}".format(bam=bam), header=False):
        if int(toks[4]) < min_map_q: continue
        # cheat and keep only perfectly matched reads (only 'M' in cigar)
        if ['M'] != [x for x in toks[5] if not x.isdigit()]: continue

        pos, flag = int(toks[3]), int(toks[1])
        ref_seq = faseq(reference, toks[2], pos,
                    pos + int("".join(x for x in toks[5] if x.isdigit())) - 1)
        aln_seq = toks[9]

        convertible = "C" if "YD:Z:f" in toks[11:] else "G"
        if convertible != "G": continue

        cpgs = [i for i, (a, b) in enumerate(zip(ref_seq, ref_seq[1:])) if a +
                b == "CG" and i < len(aln_seq)]
        # from faseq() we get 1 base on either end of seq, so here we string
        # back to the original
        ref_seq = ref_seq[1:]
        second_read = int(bool(flag & 0x40))
        for idx in cpgs:
            ref, seq = ref_seq[idx], aln_seq[idx]
            assert ref == "G", ref
            if ref == "G" and not seq in "GA": continue
            meth_counts[second_read][idx][seq] += 1

        n_to_check -= 1
        if n_to_check == 0: break

    read_length = max(meth_counts[0].keys())
    bias_list = []
    fh.write("base\tread_1\tread_2\n")
    for i in range(read_length + 1):
        ratios = []
        for read in (0, 1): # paired end
            mi = meth_counts[read][i]
            ratios.append((mi["G"]) / (float(mi["G"] + mi["A"]) or 1.))
        bias_list.append(dict(base=i + 1, read_1=ratios[0], read_2=ratios[1]))
        fh.write("%s\t%.3f\t%.3f\n" % (i + 1, ratios[0], ratios[1]))
    sys.stderr.write("wrote to %s\n" % fh.name)
    fh.close()

    try:
        plot_bias(bias_list, fpath)
    except ImportError:
        pass

def plot_bias(bias_list, fpath):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns

    sys.stderr.write("saving to %s\n" % fpath)
    p = sns.color_palette("deep", desat=.8)
    c1, c2 = p[0], p[-1]
    c1, c2 = sns.color_palette("Set1", 2)

    mpl.rc("figure", figsize=(8, 4))
    f, ax1 = plt.subplots(1)

    xs = [int(x['base']) for x in bias_list]

    ax1.plot(xs, [100 * float(x['read_1']) for x in bias_list], c=c1,
             label="read 1")
    ax1.plot(xs, [100 * float(x['read_2']) for x in bias_list], c=c2,
             label="read 2")
    ax1.axvline(x=4, c="#aaaaaa", alpha=0.8, linewidth=3, zorder=-1)
    ax1.axvline(x=max(xs) - 4, c="#aaaaaa", alpha=0.8, linewidth=3, zorder=-1)

    ax1.legend(loc='upper center')
    ax1.set_xlim(min(xs), max(xs))
    ax1.set_ylabel('mean CG % methylation')
    ax1.set_xlabel('position along read')
    ax1.set_title('Methylation Bias Plot (vertical lines at 4 bases from end)')
    ax1.grid('off')

    f.tight_layout()
    f.savefig(fpath)

if __name__ == "__main__":

    bam, ref = sys.argv[1:3]
    check_bias(bam, ref)
