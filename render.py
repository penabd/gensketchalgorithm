#!/usr/bin/python
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def case1():
    """Defines ellipse
    """
    sources = [(0,0)]
    
    # Load scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_ellipse_data.txt")

    return data, sources

def case2():
    """
    """
    sources = [(0.25,0.25), [-0.25, -0.25]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_2.txt")

    return data, sources


def case3():
    """
    """
    sources = [(1,0.5), [1, -0.5]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_3.txt")

    return data, sources


def case4():
    """
    """
    sources = [(1,1), [0,-1], [-1,1]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_4.txt")

    return data, sources


def case5():
    """
    """
    sources = [(1,1), [-2,0], [1,0], [4,0]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_6.txt")

    return data, sources

def case6():
    """
    """
    sources = [(1,1), [-2,0], [1,0], [4,0]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_7.txt")

    return data, sources

def coordinates(input):
    split = input.split(',')
    return [float(split[0][1::]),  float(split[1][:-1:])]

def main():
    testCase = "case6"
    if testCase == "case1":
        data, sources = case1()
    elif testCase == "case2":
        data, sources = case2() 
    elif testCase == "case3":
        data, sources = case3()
    elif testCase == "case4":
        data, sources = case4()
    elif testCase == "case5":
        data, sources = case5()
    elif testCase == "case6":
        data, sources = case6()
    else:
        raise ValueError("Invalid test case")
    
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:,0], 
                            data[:,1],
                            data[:,2], levels=10,
                            colors="grey",)

    color = 'k'


    with open("sketch_plot.txt", "r") as fp:
        Lines = fp.readlines()
        for line in Lines:
            split = line.split()
            if split[0] == 'Ellipse':
                [centerx, centery] = coordinates(split[1])
                width = 2 * float(split[2])
                height = 2 * float(split[3])

                ellipse = Ellipse(xy=(centerx, centery), width=width, height=height, edgecolor='forestgreen', fc='greenyellow', fill=True)
                ax.add_patch(ellipse)
            elif split[0] == 'Line':
                [startx, starty] = coordinates(split[1])
                [endx, endy] = coordinates(split[2])

                ax.plot([startx, endx], [starty, endy], color, lw=0.5)
            elif split[0] == 'Pen':
                color = split[1]

    ax.set_aspect('equal', adjustable='datalim')
    plt.colorbar(contour, ax=ax)  # Add colorbar for the contour plot
    plt.ylim(-5,5)
    pp = PdfPages('sketch.pdf')
    pp.savefig(fig)  # Save the figure with both contour and sketch
    pp.close()
    plt.show()  # Display the plot

if __name__ == '__main__':
    main()
