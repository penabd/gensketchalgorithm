#!/usr/bin/python

#from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def coordinates(input):
    split = input.split(',')
    return [float(split[0][1::]),  float(split[1][:-1:])]

def main():
    plt.figure()
    ax = plt.gca()

    color = 'k'

    with open("sketch_plot.txt", "r") as fp:
        Lines = fp.readlines()
        for line in Lines:
            split = line.split()
            if split[0] == 'Line':
                # print("Line")
                [startx, starty] = coordinates(split[1])
                [endx, endy] = coordinates(split[2])

                plt.plot([startx, endx], [starty, endy], color, lw=0.5)
            #    plt.arrow((startx + endx) / 2, (starty + endy) / 2, (endx - startx) / 10, (endy - starty) / 20, color=color, shape='full', length_includes_head=False, head_width=.01)
            elif split[0] == 'Pen':
                color = split[1]

    pp = PdfPages('sketchfinal.pdf')
    pp.savefig()
    pp.close()

if __name__ == '__main__':
    main()
