import re
def parse_time(log_file):

    for line in open(log_file):
        if "CPU time" in line:
            _, time = [x.strip() for x in line.split(":")]
            time, unit = time.split()
            assert unit == "sec."
            return float(time)
    raise Exception(log_file)

def parse_mem(log_file):

    for line in open(log_file):
        if "Max Memory" in line:
            _, mem = [x.strip() for x in line.split(":")]
            mem, unit = mem.split()
            assert unit == "MB"
            return float(mem)

    raise Exception(log_file)

def parse_prog(log_file):
    import os.path as op
    p = op.basename(log_file)
    if p.startswith('trim-'):
        p = p[5:]
    return p.split("-")[0]

def parse_sim_real(log_file):
    import os.path as op
    p = op.basename(log_file)
    if p.startswith('trim-'):
        p = p[5:]
    return p.split("-")[1].split(".")[0]

def parse_trim(log_file):
    return ["no", "yes"][int("trim" in log_file)]


FMT = "{trimmed} &    {program} & {time}(min) & {mem}(GB) & {dataset} \\\\"

def main(logs):
    fmt = FMT
    print fmt.replace("{", "").replace("}", "")
    fmt = re.sub("\(.+?\)", "", fmt)
    for log in logs:
        trimmed = parse_trim(log)
        mem = '%.2f' % (parse_mem(log) / 1000.)
        time = '%.2f' % (parse_time(log) / 60.)
        program = parse_prog(log)
        dataset = parse_sim_real(log)
        print fmt.format(**locals())

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
