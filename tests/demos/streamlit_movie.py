# A streamlit app to display a phase screen and the associated speckle pattern
import context
import numpy as np
import streamlit as st
from MegaScreen import MegaScreen
import time
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


def autoscale(screen):
    smax = np.amax(screen)
    smin = np.amin(screen)
    return (screen - smin) / (smax - smin)


def main():
    # Interactive Streamlit elements, like these sliders, return their value.
    # This gives you an extremely simple interaction model.
    iterations = st.sidebar.slider("Iterations", 10, 20_000, 1000, 100)
    Dr0 = st.sidebar.slider("D/r_0", 1.0, 20.0, 10.0, 0.5)
    diameter = st.sidebar.slider("Diameter samples", 10, 256, 100, 1)
    r0 = diameter / Dr0
    width = diameter * st.sidebar.slider("Display zoom", 1, 10, 4, 1)
    windowSize = st.sidebar.slider("Speckle window (pixels)", 10, 400, 200, 10)
    sleep = st.sidebar.slider("Sleep (milliseconds)", 0.0, 1000.0, 200.0, 1.0)
    theta = st.sidebar.slider("Wind direction (degrees from South)", 0, 360, 30, 10)

    # Non-interactive elements return a placeholder to their location
    # in the app. Here we're storing progress_bar to update it later.
    progress_bar = st.sidebar.progress(0)

    # These two elements will be filled in later, so we create a placeholder
    # for them using st.empty()
    frame_text = st.sidebar.empty()
    st.title("MegaScreen demo app")
    screen_im = st.empty()
    speckle_im = st.empty()
    st.markdown(
        """
  The upper frame shows an infinite phase screen blowing across a telescope aperture and the lower frame the resulting speckle pattern seen in the telescope focal plane when observing an unresolved target. You can change the simulation parameters in the sidebar on the left and the simulation will restart immediately with the new parameters.

  The frame rate of this movie is limited to about 5Hz (i.e. you need to set sleep >200ms) when running over the web: frames can be played at > 30Hz on a local machine.

Find out more at the [github respository](https://github.com/dbuscher/megascreen)
"""
    )
    frame_num = 0
    mask = circular_mask(diameter)
    for screen in MegaScreen(
        r0=r0,
        windowShape=(diameter, diameter),
        theta=theta * np.pi / 180,
        numIter=iterations,
    ):
        # Here were setting value for these two elements.
        frame_num += 1
        progress_bar.progress(frame_num / iterations)
        frame_text.text("Frame %i/%i" % (frame_num + 1, iterations))
        image = speckle_image(screen, mask)
        offset = (image.shape[0] - windowSize) // 2
        image = image[offset : offset + windowSize, offset : offset + windowSize]
        screen_im.image(np.where(mask, autoscale(screen), 0.5), width=width)
        speckle_im.image(autoscale(image), width=width)
        time.sleep(sleep / 1000)

    # We clear elements by calling empty on them.
    progress_bar.empty()
    frame_text.empty()

    # Streamlit widgets automatically run the script from top to bottom. Since
    # this button is not connected to any other logic, it just causes a plain
    # rerun.
    st.button("Re-run")


main()
