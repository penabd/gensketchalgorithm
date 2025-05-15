This repository contains preliminary implementations of the Sketch Algorithm. In this simulation, we don't include any physics or control laws.

Below are descriptions of what various files represent:

render.py - Renders the output from the output file (sketch_plot.txt) of one of the C++ code into a visual showing the path of the two robots. Currently overlayed on the contour plots generated from the data in the data folder. 

sketchalgorithmautogaussian.cpp - The Sketch Algorithm for Gaussians, to run this code just compile using any C++ compiler, we recommend g++ without any flags.

This should be compiled using:
clang++ -std=c++17 -g -O0 -I<path-to-eigen3> sketchalgorithmautogaussian.cpp


Details on input format:

The variable num in line 1417 indicates the number of gaussians to be included in the plume. Then gaussianVar and gaussianCenter are two vectors that
will contain the variances.

The polygon is given as a sequence of successive vertices in the file sketchinput.txt and then it outputs the result (the path of the two robots) in sketch_plot.txt 


Other Render Files:
renderfinal.py  and renderinput.py have not been modified for the general
sketch algorithm.

renderinput.py - renders a good image of a polygon and outputs the set of its vertices in successive sequence. Copy this in a file named "sketchinput.txt" after

For debugging purposes:
render.py renders gaussians to visualize the input


Summary:

To run gaussian test cases - compile and run sketchalgorithmautogaussian.cpp and then feed the output to render.py

Finally, one can tweak with the epsilon values in the main function however they wish upto the algorithm specifications (see the paper).
