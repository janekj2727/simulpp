#!/usr/bin/python3

# plot with MACSIMUS syntax using pandas plot

import matplotlib
import pandas as pd
import matplotlib.pyplot as plot
import matplotlib.ticker as ticker
import sys
import numpy as np

# plot.rcParams["mathtext.fontset"] = "stixsans"
# plot.rcParams["text.usetex"] = True


def print_help():
    print('Plot using MACSIMUS syntax and pandas.plot()')
    print('For data syntax, check MACSIMUS plot')
    print('Usage:')
    print('    plot.py [OPTIONS] DATA')
    print('Where options are (% are replaced by spaces):')
    print('-b ... block data (e.g. -b100 = block by 100)')
    print('-l ... legend (e.g. -l["Etot","Ekin"])')
    print('-p ... legend position (e.g. -p["lower%left"])')
    print('-t ... title (e.g. -t["Graph%title"])')
    print('-x ... xaxis label (e.g. -x["time(s)"])')
    print('-y ... yaxis label (e.g. -y["energy(J)"])')
    print('-g ... grid off (default grid on)')
    print('-e ... save figure as .pdf (e.g. -e[plot3.pdf])')
    print('-f ... set global font size (e.g. -f12), default 10')
    print('-s ... set plot size (default -s[7,4])
    quit()


def parse_option(option):
    if (option[1] == 'h'):
        print_help()
    elif (option[1] == 'b'):
        global block
        block = int(option[2:])
    elif (option[1] == 'g'):
        global grd
        grd = False
    elif (option[1] == 'l'):
        global leg
        leg = option[2:].strip('[]').replace('%',' ').split(',')
    elif (option[1] == 't'):
        global tit
        tit = str(option[2:].strip('[]')).replace('%',' ')
    elif (option[1] == 'p'):
        global lloc
        lloc = str(option[2:].strip('[]')).replace('%',' ')
    elif (option[1] == 'x'):
        global xlab
        xlab = str(option[2:].strip('[]')).replace('%',' ')
        print("xlabel = " + xlab)
    elif (option[1] == 'y'):
        global ylab
        ylab = str(option[2:].strip('[]')).replace('%',' ')
    elif (option[1] == 's'):
        global siz
        siz = list(map(int,option[2:].strip('[]').split(',')))
    elif (option[1] == 'e'):
        global savename
        savename = str(option[2:].strip('[]'))
    elif (option[1] == 'f'):
        global fontsize
        fontsize = int(option[2:])
    else:
        print("Unknown option")


def parse_string(string, currfile):
    for c in range(65, 85):
        string = string.replace(chr(c), "df[" + str(currfile) + "][" + str(c - 65) + "]")
    return string


if (len(sys.argv) == 1) or (sys.argv[1] == '-h'):
    print_help()

print(sys.argv)

fileopened = False
df = []
currfile = int(-1)
currcolumn = 0
dffinal = pd.DataFrame()
block = 1
leg = []
lloc = None
grd = True
tit = None
xlab = None
ylab = None
siz = [7,4]
savename = None
fontsize = 0
labfontsize = 14


for arg in sys.argv[1:]:
    if (arg[0] == '-'):
        parse_option(arg)
    else:
        argsplit = arg.split(':')
        if (len(argsplit) > 2):
            indexcol = None
            if (int(argsplit[1]) > 0):
                indexcol = int(argsplit[1]) - 1
            print("Opening file: " + argsplit[0])
            df.append(pd.read_csv(
                str(argsplit[0]), sep="\s+", comment="#", header=None, index_col=indexcol,error_bad_lines=False))
            currfile += 1
            fileopened = True
            if (not argsplit[2].isdigit()):
                colstr = parse_string(argsplit[2].strip('"'), currfile)
                # print("here")
            else:
                colstr = "df[" + str(currfile) + "][int(argsplit[2])-1]"
            dffinal[currcolumn] = eval(colstr)
            currcolumn += 1
        elif (fileopened == True):
            if (not argsplit[-1].isdigit()):
                colstr = parse_string(argsplit[-1], currfile)
            else:
                colstr = "df[" + str(currfile) + "][int(argsplit[2])-1]"
            dffinal[currcolumn] = eval(colstr)
            currcolumn += 1
        else:
            print("ERROR")

print(dffinal)

if (fontsize > 0):
    matplotlib.rcParams.update({'font.size':fontsize})
    labfontsize=1.5*fontsize

dff = dffinal.rolling(block).mean()
dff.plot()
plot.grid(grd)


ax = plot.gca()
ax.ticklabel_format(useOffset=False)
if (lloc != None):
    ax.legend(leg,loc=lloc)
else:
    ax.legend(leg)
ax.set_title(tit)
ax.set_xlabel(xlab,usetex=True,fontsize=labfontsize)
ax.set_ylabel(ylab,usetex=True,fontsize=labfontsize)

fig = plot.gcf()
fig.set_size_inches(siz)
fig.tight_layout()


plot.show()
print(dff)

if (savename != None):
    fig.savefig(savename,format='pdf')
