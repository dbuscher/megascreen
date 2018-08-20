# Make movie of phase screen

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import MegaScreen
import functools


def test1(diameter=128, screenSize=256, numIter=100, dx=3.3, dy=10):
    tileGenerator=MegaScreen.SplineTiles(MegaScreen.McGlameryScreen(r0=10, L0=100,
                                                                    nfft=256))
    screenGenerator=MegaScreen.SlidingWindows(tileGenerator,(diameter,diameter),dx=3.3,theta=np.pi/3)
    ScreenMovie(screenGenerator)

def test(r0=10,L0=1e4,diameter=128,dx=3.3):
    generator=MegaScreen.MegaScreen(r0,L0,windowShape=[diameter,diameter],dx=dx,theta=np.pi/3)
    ScreenMovie(generator)

def test3(r0=10,L0=1e4,diameter=128,dx=3.3):
    ScreenMovie(SingleGenerator(r0,L0,diameter,dx,which=0))

def SingleGenerator(r0,L0,diameter,dx,which=0):
    for screens in MegaScreen.MegaScreen(r0,L0,windowShape=[diameter,diameter],dx=dx,debug=True):
        yield screens[which]

def ScreenMovie(screenGenerator):
    fig = plt.figure()
    screen = next(screenGenerator)
    im = plt.imshow(screen, cmap=plt.cm.gray,interpolation='none',animated=True)
    ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True, fargs=(im,screenGenerator))
    plt.show()

def updatefig(i,im,screenGenerator):
    screen=next(screenGenerator)
    im.set_data(screen)
    im.autoscale()
    return im,


# Generic command-line interface to run test script or give full control
if __name__ == '__main__':
    if len(sys.argv) == 1:
        test()
    else:
        exec(sys.argv[1])
