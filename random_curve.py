# File name: random_curve.py
# Author: Marcos Antonio Avila
# Date created: 10/23/2019
# Date last modified: 10/24/2019
# Python Version: 3.7

import matplotlib.pyplot as plt
import numpy as np

# Image file name
FILE_NAME = "curby.png"
# Number of points the curve will touch
VIA_POINT_COUNT = 5
# Number of random points to draw on curve
RAND_POINT_COUNT = 2
# Time required to reach next point
STEP_SIZE = 10

# Bounds for x and y positions (via points)
MIN_X_POS = 0
MAX_X_POS = 800
MIN_Y_POS = 0
MAX_Y_POS = 600
# Bounds for instantaneous x and y velocities at each via point
VEL_THRESH = 500
# Bounds for instantaneous x and y accelerations at each via point
ACC_THRESH = 100

# Randomly generate x and y positions
xPos = np.random.randint(MIN_X_POS, MAX_X_POS, VIA_POINT_COUNT)
yPos = np.random.randint(MIN_Y_POS, MAX_Y_POS, VIA_POINT_COUNT)
# Randomly generate x and y velocities
xVel = np.random.randint(-VEL_THRESH, VEL_THRESH, VIA_POINT_COUNT)
yVel = np.random.randint(-VEL_THRESH, VEL_THRESH, VIA_POINT_COUNT)
# Randomly generate x and y accelerations
xAcc = np.random.randint(-ACC_THRESH, ACC_THRESH, VIA_POINT_COUNT)
yAcc = np.random.randint(-ACC_THRESH, ACC_THRESH, VIA_POINT_COUNT)

# Create a matrix of coefficients for 5th order functions which connect continuously at via points
xCoeffMat = [tuple([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])] * VIA_POINT_COUNT
yCoeffMat = [tuple([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])] * VIA_POINT_COUNT

# Coefficient calculation for x and y parametric functions
for index in range(VIA_POINT_COUNT):
    xCoeffMat[index] = (xPos[index],
                        xVel[index],
                        xAcc[index] / 2,
                        ((20*(xPos[(index + 1) % VIA_POINT_COUNT] - xPos[index]) -
                         STEP_SIZE * (8 * xVel[(index + 1) % VIA_POINT_COUNT] + 12 * xVel[index]) -
                         STEP_SIZE ** 2 * (3 * xAcc[index] - xAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 3)),
                        ((30 * (xPos[index] - xPos[(index + 1) % VIA_POINT_COUNT]) +
                         STEP_SIZE * (14 * xVel[(index + 1) % VIA_POINT_COUNT] + 16 * xVel[index]) +
                         STEP_SIZE ** 2 * (3 * xAcc[index] - 2 * xAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 4)),
                        ((12 * (xPos[(index + 1) % VIA_POINT_COUNT] - xPos[index]) -
                         STEP_SIZE * 6 * (xVel[(index + 1) % VIA_POINT_COUNT] + xVel[index]) -
                         STEP_SIZE ** 2 * (xAcc[index] - xAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 5)))
    yCoeffMat[index] = (yPos[index],
                        yVel[index],
                        yAcc[index] / 2,
                        ((20*(yPos[(index + 1) % VIA_POINT_COUNT] - yPos[index]) -
                         STEP_SIZE * (8 * yVel[(index + 1) % VIA_POINT_COUNT] + 12 * yVel[index]) -
                         STEP_SIZE ** 2 * (3 * yAcc[index] - yAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 3)),
                        ((30 * (yPos[index] - yPos[(index + 1) % VIA_POINT_COUNT]) +
                         STEP_SIZE * (14 * yVel[(index + 1) % VIA_POINT_COUNT] + 16 * yVel[index]) +
                         STEP_SIZE ** 2 * (3 * yAcc[index] - 2 * yAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 4)),
                        ((12 * (yPos[(index + 1) % VIA_POINT_COUNT] - yPos[index]) -
                         STEP_SIZE * 6 * (yVel[(index + 1) % VIA_POINT_COUNT] + yVel[index]) -
                         STEP_SIZE ** 2 * (yAcc[index] - yAcc[(index + 1) % VIA_POINT_COUNT])) /
                        (2 * STEP_SIZE ** 5)))

# Create a list of time slices from 0 until time it takes to reconnect curve (VIA_POINT_COUNT * STEP_SIZE)
t = np.linspace(0, VIA_POINT_COUNT * STEP_SIZE, int(VIA_POINT_COUNT * STEP_SIZE * 100))

# Create lists to hold intermediate x and y positions
x = [0.0] * len(t)
y = [0.0] * len(t)

# Create list to hold euclidean distances between adjacent intermediate coordinates
t_dists = [0.0] * len(t - 1)
# Create variable to hold approximated total distance
dist = 0.0

# Iterate through all time slices
for t_index in range(len(t)):
    # Get spline section index
    s_index = int(t[t_index] // STEP_SIZE)

    # Get coefficients of appropriate spline section
    xCoeff = xCoeffMat[s_index % VIA_POINT_COUNT]
    yCoeff = yCoeffMat[s_index % VIA_POINT_COUNT]

    # Calculate x and y intermediate values
    x[t_index] = (xCoeff[0] +
                  xCoeff[1] * (t[t_index] - (s_index * STEP_SIZE)) +
                  xCoeff[2] * (t[t_index] - (s_index * STEP_SIZE))**2 +
                  xCoeff[3] * (t[t_index] - (s_index * STEP_SIZE))**3 +
                  xCoeff[4] * (t[t_index] - (s_index * STEP_SIZE))**4 +
                  xCoeff[5] * (t[t_index] - (s_index * STEP_SIZE))**5)
    y[t_index] = (yCoeff[0] +
                  yCoeff[1] * (t[t_index] - (s_index * STEP_SIZE)) +
                  yCoeff[2] * (t[t_index] - (s_index * STEP_SIZE))**2 +
                  yCoeff[3] * (t[t_index] - (s_index * STEP_SIZE))**3 +
                  yCoeff[4] * (t[t_index] - (s_index * STEP_SIZE))**4 +
                  yCoeff[5] * (t[t_index] - (s_index * STEP_SIZE))**5)

    # Skip first distance calculation, since a previous value has not been generated
    if t_index == 0:
        continue

    # Calculate euclidean distance from one intermediate value to the next
    t_dists[t_index - 1] = np.sqrt((x[t_index] - x[t_index - 1])**2 + (y[t_index] - y[t_index - 1])**2)
    # Increment total distance by previously calculated intermediate distance
    dist += t_dists[t_index - 1]

# Generate random points on curve
for p_index in range(RAND_POINT_COUNT):
    # Get random distance from 0 to total distance
    r_dist = dist * np.random.rand()

    # Iterate through t_dists until you reach the random distance
    for t_index in range(len(t_dists)):
        # Decrement the random distance by intermediate distance
        r_dist -= t_dists[t_index]

        # Once random distance is reached, plot x and y coordinate
        if r_dist < 0:
            plt.plot(x[t_index], y[t_index], 'ro')
            break

# Display via points
# plt.plot(xPos, yPos, 'bo')

# Display plot of x and y
plt.plot(x, y, 'k-')
plt.axis('off')
plt.savefig(FILE_NAME, bbox_inches='tight')
plt.show()
