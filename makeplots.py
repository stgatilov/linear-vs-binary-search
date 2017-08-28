import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys, os, glob, re, math


results = {}
def AddRes(name, mem, N, eltime):
    global results
    results.setdefault(mem, {}).setdefault(name, {})[N] = eltime


def ReadRes(fn):
    with open(fn, "rt") as f:
        data = f.read()

    match = re.search(r"Arrays: (\d*) x (\d*)", data)
    N = int(match.group(2))
    match = re.search(r"Memory: (\d*)", data)
    mem = int(match.group(1))

    for match in re.finditer(r"\s*([0-9.]+)\s*ns\s*:\s*(\S+)", data):
        eltime = float(match.group(1))
        name = match.group(2)
        AddRes(name, mem, N, eltime)


for fn in glob.glob("res/*.log"):
    ReadRes(fn)

# plt.loglog([1,2,3,4], [1,4,9,16], 'bo', [1,2,3,4], [16,9,9,10], 'ro', basex=2, basey=2, linestyle='-')
# plt.show()

styles = ['yx', 'rx', 'r+', 'mx', 'm+', 'k.', 'ko', 'bo', 'bs', 'yo', 'g*', 'gP', 'gd', 'm*', 'c*']

dpi = 150

for mem, graphs in results.items():
    args = []
    names = []
    argsPE = []
    argsLog = []

    idx = 0
    for name, graph in graphs.items():
        if ('linear' in name and 'scalar' in name):
            continue
        X = []
        Y = []
        Z = []
        W = []
        for N, eltime in graph.items():
            X.append(N)
            Y.append(eltime)
            Z.append(eltime / N)
            W.append(eltime / math.log(N, 2.0))
        args += [X, Y, styles[idx]]
        argsPE += [X, Z, styles[idx]]
        argsLog += [X, W, styles[idx]]
        names.append(name)
        idx += 1
        print("%s: %s" % (name, args[-1]))


    title = "(memory = %dB)" % mem
    if len(sys.argv) > 1:
        title = sys.argv[1] + " " + title

    ax = plt.axes()
    ax.set_title(title)
    ax.loglog(*args, basex=2, basey=2, linestyle='-')
    ax.set_xlabel("Array length (N)")
    ax.set_ylabel("Time per search, ns")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", color='0.8', linestyle=':')
    ax.legend(names, loc=2, prop={'size': 6})
    ax.get_yaxis().get_minor_locator().subs([1.25, 1.5, 1.75])
    ax.get_yaxis().set_minor_formatter(ticker.FuncFormatter(lambda x,p: str(int(x))))
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    #plt.show()
    plt.savefig('res/plot_search_%d.png' % mem, bbox_inches='tight', dpi=dpi)
    plt.gcf().clear()

    ax = plt.axes()
    ax.set_title(title)
    ax.semilogx(*argsPE, basex=2, linestyle='-')
    ax.set_xlabel("Array length (N)")
    ax.set_ylabel("Time per element, ns")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", color='0.8', linestyle=':')
    ax.legend(names, loc=1, prop={'size': 6})
    ax.set_ylim(0.0, 0.5)
    ax.get_yaxis().set_minor_locator(ticker.MultipleLocator(0.01))
    ax.get_yaxis().tick_right()
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    #plt.show()
    plt.savefig('res/plot_elem_%d.png' % mem, bbox_inches='tight', dpi=dpi)
    plt.gcf().clear()

    ax = plt.axes()
    ax.set_title(title)
    ax.semilogx(*argsLog, basex=2, linestyle='-')
    ax.set_xlabel("Array length (N)")
    ax.set_ylabel("Time per one bin.search comparison, ns")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", color='0.8', linestyle=':')
    ax.legend(names, loc=2, prop={'size': 6})
    ax.set_ylim(1.0, 7.0)
    ax.get_yaxis().set_minor_locator(ticker.MultipleLocator(0.5))
    ax.get_yaxis().tick_right()
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    #plt.show()
    plt.savefig('res/plot_log_%d.png' % mem, bbox_inches='tight', dpi=dpi)
    plt.gcf().clear()
