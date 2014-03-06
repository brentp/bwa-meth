"""
runner for gsnap on methylation data.
gmap_build, gsnap and samtools must be on the PATH of the calling environment
"""

import argparse
import sys
import os
import os.path as op
import subprocess as sp

def sh(cmd, log=sys.stderr, wait=True):
    print >>sys.stderr, "[running command] %s" % cmd
    p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    if wait:
        for line in p.stderr:
            print >>log, line,
        p.wait()
        if p.returncode != 0 or "aborted" in p.stderr.read():
            sys.exit(p.returncode)
    return p

def nopen(f, mode="rb"):
    if not isinstance(f, basestring): return f
    return open(f, mode)


def gmap_built(ref_dir, ref_base):
    for ext in ("chromosome", "contig.iit", "ref12153offsetscomp",
                "chromosome.iit", "genomecomp",
                "chrsubset", "maps", "version", "contig",
                "ref12153gammaptrs"):

        f = "%s/%s/%s.%s" % (ref_dir, ref_base, ref_base, ext)
        if not op.exists(f):
            print >>sys.stderr, "%s does not exist" % f
            return False
    return True

def cmetindexed(ref_dir, ref_base, kmer):

    for ext in ("metct12%i3offsetscomp",
                "metct12%i3gammaptrs", "metga12%i3offsetscomp",
                "metga12%i3gammaptrs"):
        ext = ext % kmer
        f = "%s/%s/%s.%s" % (ref_dir, ref_base, ref_base, ext)
        if not op.exists(f):
            print >>sys.stderr, "%s does not exist" % f
            return False
    return True

def gsnap_index(reference, kmer=15):
    ref_dir, ref_base = check_reference(reference)
    ref_base += "." + str(kmer)

    cmd = "gmap_build --no-sarray -k %(kmer)i -D %(ref_dir)s -d %(ref_base)s %(reference)s"
    cmd %= locals()
    if not gmap_built(ref_dir, ref_base):
        sh(cmd)
    else:
        print >>sys.stderr, "[skipping command] %s" % cmd

    cmd_cmet = "cmetindex -d %(ref_base)s -F %(ref_dir)s -k %(kmer)d\n"
    cmd_cmet %= locals()
    if not cmetindexed(ref_dir, ref_base, kmer):
        sh(cmd_cmet)
    else:
        print >>sys.stderr, "[skipping command] %s" % cmd

def check_reference(reference):
    ref_dir = op.dirname(reference)
    ref_base = op.splitext(op.basename(reference))[0] # locals
    assert os.access(reference, os.R_OK), ("reference not found / readable")
    assert os.access(ref_dir, os.W_OK), ("%s must be writable" % (ref_dir))
    return ref_dir, ref_base

def gsnap_meth(reference, reads, prefix, kmer=15, stranded=False,
        extra_args="", threads=1, rg=None):

    ref_dir, ref_base = check_reference(reference)
    ref_base += "." + str(kmer)
    mode = ["cmet-nonstranded", "cmet-stranded"][int(stranded)]
    reads_str = " ".join(reads)
    if any(r.endswith(".gz") for r in reads):
        extra_args = extra_args + " --gunzip"
    cmd_gsnap = "set -eo pipefail;"
    cmd_gsnap += "gsnap -B 4 --npaths 1 --quiet-if-excessive \
        --nthreads {threads} \
        -A sam -k {kmer} -D {ref_dir} -d {ref_base} --mode {mode} \
        --use-sarray 0 \
        --pairexpect 300 \
        --pairdev 250 \
        --read-group-id {rg} --read-group-name {rg} \
         {extra_args}  {reads_str}"
    cmd_gsnap += "| samtools view -bS - | samtools sort - {prefix}"
    cmd_gsnap = cmd_gsnap.format(**locals())
    sh(cmd_gsnap)

    cmd_index = "samtools index {prefix}.bam".format(**locals())
    sh(cmd_index)

def run(args):
    if args.threads is None:
        import multiprocessing
        threads = multiprocessing.cpu_count() # locals
    else:
        threads = args.threads
    reference, kmer = op.abspath(args.reference), args.kmer
    reads = map(op.abspath, args.reads)
    rg = rname(reads[0], reads[1])
    gsnap_meth(reference, reads, args.prefix, args.kmer, args.stranded,
               args.extra_args, threads, rg=rg)

def rname(fq1, fq2):
    def name(f):
        n = op.basename(op.splitext(f)[0])
        if n.endswith('.fastq'): n = n[:-6]
        if n.endswith(('.fq', '.r1', '.r2')): n = n[:-3]
        return n
    return "".join(a for a, b in zip(name(fq1), name(fq2)) if a == b) or 'bm'


def main():
    if len(sys.argv) > 1 and "index" == sys.argv[1]:
        sys.exit(gsnap_index(sys.argv[2], kmer=int(sys.argv[3])))

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", help="reference fasta")
    p.add_argument("-k", dest="kmer", help="kmer length for gsnap. default:" \
                   + "%(default)i", type=int, default=15)
    p.add_argument("-t", "--threads", dest="threads", default=None,
            help="number of threads for gsnap to use", type=int)
    p.add_argument("--stranded", action="store_true", default=False, help=\
                   "by default, non-stranded library prep is assumed")
    p.add_argument("--prefix", help="prefix for output bam",
                  default="gsnap-meth")
    p.add_argument("--extra-args", dest="extra_args", help="extra arguments"
                " to send to gsnap", default="")

    p.add_argument('reads', nargs='+', help='reads files (if 2 files are'
                   'they are assumed to be paired end)')

    try:
        args = p.parse_args()
    except:
        p.print_help()
        raise

    if (not len(args.reads) in (1, 2)) or args.reference is None:
        sys.exit(not p.print_help())

    run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
