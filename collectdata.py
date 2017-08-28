import re, os, glob

def CollectData(N, mem):
    for fn in glob.glob("search_tmp.*"):
        os.remove(fn)

    with open("search.cpp", "rt") as f:
        src = f.read()
    ints = mem / 8
    src = re.sub(r"const int SIZE = (\d*);", r"const int SIZE = %d;" % N, src)
    src = re.sub(r"const int ARR_SAMPLES = (\(.*?\))", r"const int ARR_SAMPLES = %d" % ints, src)
    src = re.sub(r"const int KEY_SAMPLES = (\(.*?\))", r"const int KEY_SAMPLES = %d" % ints, src)
    with open("search_tmp.cpp", "wt") as f:
        f.write(src)

    with open("c.bat", "rt") as f:
        bat = f.read()
    bat = bat.replace("search.cpp", "search_tmp.cpp")
    os.system(bat)

    logname = "res_%04d_%d.log" % (N, mem)
    os.system("search_tmp >res/" + logname)
    os.system("search_tmp >res/" + logname)


for fn in glob.glob("res/*"):
    os.remove(fn)

sizes = [16, 32, 64, 128, 256, 512, 1024]
#sizes = [128, 256, 512, 1024, 2048, 4096]
for s in sizes:
    CollectData(s, 64<<10)
#    CollectData(s, 512<<10)
