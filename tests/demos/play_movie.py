#!/usr/bin/env python3
# Make movie of phase screen and show in a local window
import context
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
from MegaScreen import MegaScreen
from matplotlib.animation import FuncAnimation


def screen_movie(
    r0=10, L0=1000, diameter=200, dx=3.5, theta=np.radians(20.0), num_sigma=2
):
    limit = num_sigma * np.sqrt(0.0863) * (L0 / r0) ** (5 / 6)
    im = plt.imshow(
        np.zeros((diameter, diameter)),
        cmap=plt.cm.gray,
        interpolation="none",
        animated=True,
        vmin=-limit,
        vmax=limit,
    )
    yield im,  # Dummy first frame to draw axes etc
    for screen in MegaScreen(
        r0=r0, L0=L0, windowShape=(diameter, diameter), dx=dx, theta=theta
    ):
        im.set_data(screen)
        yield im,


if __name__ == "__main__":
    fig = plt.figure()
    anim = FuncAnimation(
        fig, lambda x: x, frames=screen_movie(), interval=30, blit=True
    )
    plt.show()
