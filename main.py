import math
import numpy as np
import plotly.graph_objects as go   # plotly is much faster at 3D scatter than matplotlib
import random
import scipy.special as sci         # to swiftly compute the wave function
import time

# TODO: Program an GUI, probably using PyQT.
print("========== HYDROGEN-LIKE ATOM 3D PROBABILITY DENSITY GENERATOR ==========")
print("                         created by Riley Reese                          \n")

principal_qn = int(input("Generate probability densities up to what principal quantum number 'n'? "))
maximum = int(input("What is the maximum x, y, and z of the bounding box? "))
minimum = int(input("What is the minimum x, y, and z of the bounding box? "))
scale = float(input("What is the scale? "))

# TODO: This is messy. Refactor - if I can.
# Define the principle, angular momentum, and magnetic quantum numbers as ranges
angular_momentum_qn = principal_qn - 1
magnetic_qn = range(-angular_momentum_qn, angular_momentum_qn+1)
angular_momentum_qn = range(angular_momentum_qn+1)
principal_qn = range(1, principal_qn+1)

# psi, the wave function equation, the complex square root of the probability density
def wfc_equation(n, l, m, r, theta, phi):

    rho_factor = (2 * r) / (n) # named after the greek letter which represents it
    condon_shortley = (-1) ** m # this is convention
    coeffecient = np.e ** (-(rho_factor / 2)) * rho_factor ** l

    # NOTE: Originally, I had actually programmed a calculator for the associated laguerre polynomial
    # and the spherical harmonic by hand. Now that I understand both of them, however, I have swapped
    # out my painfully inefficient algorithms with the much faster scipy equivelants.

    laguerre = sci.assoc_laguerre(rho_factor, n - l - 1, 2 * l + 1,) # the associated laguerre polynomial (capital L)
    harmonic = sci.sph_harm(m, l, theta, phi) # the spherical harmonic (capital Y)

    # I actually have no idea why you multiply by root 2, but some other guy did it, so...
    return condon_shortley * coeffecient * laguerre * harmonic * (2 ** (1/2))


# Convert the real component of the complex probability from 3D polar coordinates to cartesian
def cartesian_probability_real(n, l, m, x, y, z):

    r = (x**2 + y**2 + z**2) ** (1/2)
    theta = np.arctan2(y, x)
    phi = np.arctan2((x**2 + y**2) ** (1/2), z)

    return np.absolute(np.real(wfc_equation(n, l, m, r, theta, phi)))**2


# Convert the imaginary component of the complex probability from 3D polar coordinates to cartesian
def cartesian_probability_imaginary(n, l, m, x, y, z):

    r = (x**2 + y**2 + z**2) ** (1/2)
    theta = np.arctan2(y, x)
    phi = np.arctan2((x**2 + y**2) ** (1/2), z)

    return np.absolute(np.imag(wfc_equation(n, l, m, r, theta, phi)))**2


# Populate a set of uniform points, which will eventually be plotted on to a
# 3D scatter graph. At each point, calculate the probability density.
def create_points(mini, maxi, s, n, l, m):

    points_x = []
    points_y = []
    points_z = []
    probabilities = []
    colors = []
    i = 0

    # List through every integer point in space within the bounding box
    for x in range(mini, maxi):
        for y in range(mini, maxi):
            for z in range(mini, maxi):
                # Generate that point's electron probability
                p = cartesian_probability_real(
                    n, l, m,
                    x * s, y * s, z * s)

                # Add that point to the list
                # TODO: 4D Array instead of 4 Lists
                probabilities.append(p)
                points_x.append(x * s)
                points_y.append(y * s)
                points_z.append(z * s)

    # calculate the maximum of p for normalization, exlcuding 0 (such that no div by 0 errors)
    maxp = max(p for p in probabilities if p != 0)
    meanp = sum(probabilities) / len(probabilities)

    # normalize all p values to the interval [0, 1]
    for p in range(len(probabilities)):
        probabilities[p] = probabilities[p] / maxp

    # generate color pallete
    # TODO: Improve visualization of the probability distribution.
    for p in range(len(probabilities)):
        colors.append((probabilities[p], 0, 1 - probabilities[p], 1 / 3 * math.sqrt(probabilities[p])))

    # The function generates one one Scatter3D graph for its visualization
    # TODO: Optional generation of a Scatter2D visualization at z = 0 on the positive x/y plane.
    return go.Scatter3d(visible = False,
                        name = f"{n_qn} {l_qn} {m_qn}",
                        mode='markers',
                        x = points_x,
                        y = points_y,
                        z = points_z,
                        marker = dict(color=colors))

# The marker data list holds all of the Scatter3D graphs
marker_data = []
total_time = 0
configurations = 0

# Generate the probability distribution for each possible configuration at n
# TODO: The lqn and mqn if statements are inelegant. Find another way.
# TODO: The probability distribution generator should be a function
# TODO: Many probability distributions are just reorientations. The potential complexity reduction is dramatic if I apply that.

print("Generating probability distributions...\n")
for n_qn in principal_qn:
    for l_qn in (range(0,n_qn) if range(0,n_qn) else range(1)):
        for m_qn in (range(-l_qn, l_qn+1) if range(-l_qn, l_qn+1) else range(1)):

            print(f"\tGenerating probability distribution {n_qn} {l_qn} {m_qn}...")
            t_initial = time.time()
            marker_data.append(create_points(minimum, maximum, scale, n_qn, l_qn, m_qn))
            t_final = time.time()
            delta_t = t_final - t_initial
            total_time += delta_t
            configurations += 1
            print(f"\tGenerated probability distribution {n_qn} {l_qn} {m_qn} in {round(delta_t, 2)} seconds!\n")

print(f"Generated {configurations} configurations in {round(total_time, 1)} seconds.")
print("\nLoading scatter plots...")
graph_t1 = time.time()

marker_data[0].visible = True
fig = go.Figure()
fig.add_traces(marker_data)

# weird slider code
# TODO: Refactor the slider code
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "Configuration: " + str(i)}],  # layout attribute
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue={"prefix": "Configuration: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    sliders=sliders
)


fig.show()
graph_t2 = time.time()
graph_deltat = round(graph_t2 - graph_t1, 1)
print(f"Succesfully loaded scatter plots in {graph_deltat} seconds!\n")
