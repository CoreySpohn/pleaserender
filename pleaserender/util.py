import astropy.units as u
import numpy as np


def filter_data(dataset, animation_value, animation_key, animation_style):
    match animation_style:
        case "Cumulative":
            return dataset.where(dataset[animation_key] <= animation_value, drop=True)
        case "Single point":
            return dataset.where(dataset[animation_key] == animation_value, drop=True)
        case "Trailing":
            # Data of trailing 10 points
            # Finding the index of the specified animation value
            current_index = np.argmax(dataset[animation_key].values == animation_value)
            # current_index = dataset[animation_key].argmax()
            # min_index = 0  # Assuming 0 is the starting index in your dimension
            return dataset.isel({animation_key:slice(max(0, current_index-9), current_index+1)})

            # if current_index - 10 > min_index:
            #     return dataset.where(
            #             (dataset[animation_key].values <= current_index) &
            #             (dataset[animation_key].values >= current_index - 10),
            #             drop=True
            #             )
            # else:
            #     breakpoint()
            #     return dataset.where(
            #             dataset[animation_key].values <= current_index, drop=True
            #             )

def calc_observer_position(elev, azim, dist):
    """
    This is a simple function to calculate the observer's position in Cartesian
    coordinates given the observer's elevation and azimuth angles and distance
    from the origin. I think it is correct but re-use at your own risk.
    """
    # Convert angles to radians
    elev_rad = np.deg2rad(elev)
    azim_rad = np.deg2rad(azim)

    x = dist*np.cos(elev_rad) * np.cos(azim_rad)
    y = dist*np.cos(elev_rad) * np.sin(azim_rad)
    z = dist*np.sin(elev)

    return u.Quantity([x, y, z])

def lims_and_ticks(ax, axes, val0, valf, dval, offset, use_minor=True):
    if 'x' in axes:
        ax.set_xlim([val0-offset, valf+offset])
        ax.set_xticks(np.arange(val0, valf+dval/4, dval))
        if use_minor:
            ax.set_xticks(np.arange(val0+dval/2, valf+dval/4, dval), minor=True)
    if 'y' in axes:
        ax.set_ylim([val0-offset, valf+offset])
        ax.set_yticks(np.arange(val0, valf+dval/4, dval))
        if use_minor:
            ax.set_yticks(np.arange(val0+dval/2, valf+dval/4, dval), minor=True)
    if 'z' in axes:
        ax.set_zlim([val0-offset, valf+offset])
        ax.set_zticks(np.arange(val0, valf+dval/4, dval))
        if use_minor:
            ax.set_zticks(np.arange(val0+dval/2, valf+dval/4, dval), minor=True)
    return ax
