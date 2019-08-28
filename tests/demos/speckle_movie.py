#!/usr/bin/env python3
# Display an endless movie of the evolution of a simulated speckle pattern
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
from MegaScreen import MegaScreen
from numpy.fft import fft2, fftshift


def circular_mask(gridSize, diameter=None):
    if diameter == None:
        diameter = gridSize
    radius = diameter / 2.0
    y, x = (np.mgrid[0:gridSize, 0:gridSize] - (gridSize - 1.0) / 2.0) / radius
    temp = x + 1j * y
    r = np.abs(temp)
    return np.less_equal(r, 1.0)


def speckle_image(screen, mask, oversample=4):
    pupil = np.exp(1j * screen) * mask
    image = fftshift(abs(fft2(pupil, s=oversample * np.array(screen.shape)))) ** 2
    return image


def speckle_movie(diameter=100, windowSize=50, **kwargs):
    mask = circular_mask(diameter)
    im = plt.imshow(
        np.zeros((windowSize, windowSize)),
        cmap=plt.cm.gray,
        interpolation="none",
        animated=True,
    )
    yield im,  # Dummy first frame to draw axes etc
    for screen in MegaScreen(windowShape=(diameter, diameter), **kwargs):
        image = speckle_image(screen, mask)
        offset = (image.shape[0] - windowSize) // 2
        im.set_data(image[offset : offset + windowSize, offset : offset + windowSize])
        im.autoscale()
        yield im,


if __name__ == "__main__":
    fig = plt.figure()
    anim = FuncAnimation(
        fig, lambda x: x, frames=speckle_movie(windowSize=200), interval=30, blit=True
    )
    plt.show()
