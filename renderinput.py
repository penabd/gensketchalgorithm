import cv2
import numpy as np

# load image
img = cv2.imread('polygonbasic2.png')

gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# apply canny edge detection
edges = cv2.Canny(gray, 60, 160)

# apply morphology close
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
morph = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)

# get contours and keep largest
contours = cv2.findContours(morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
contours = contours[0] if len(contours) == 2 else contours[1]
big_contour = max(contours, key=cv2.contourArea)

  # draw contour
contour = img.copy()
cv2.drawContours(contour, [big_contour], 0, (0,0,255), 1)

# get number of vertices (sides)
peri = cv2.arcLength(big_contour, True)
approx = cv2.approxPolyDP(big_contour, 0.01 * peri, True)
print('number of sides:',len(approx))

# get the sides of the polygon
for vertex in approx:
    print (vertex)

