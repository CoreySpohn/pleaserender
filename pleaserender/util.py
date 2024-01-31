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

def calc_viewer_position(elev, azim, dist):
    """
    This is a simple function to calculate the viewer's position in Cartesian
    coordinates given the viewer's elevation and azimuth angles and distance
    from the origin. I think it is correct but re-use at your own risk.
    """
    # Convert angles to radians
    elev_rad = np.deg2rad(elev)
    azim_rad = np.deg2rad(azim)

    x = dist*np.cos(elev_rad) * np.cos(azim_rad)
    y = dist*np.cos(elev_rad) * np.sin(azim_rad)
    z = dist*np.sin(elev)

    return u.Quantity([x, y, z])

def get_nice_number(value, round=False):
    # Exponent of base 10
    exponent = np.floor(np.log10(value))

    # Fractional part
    fraction = value / 10**exponent

    if round:
        if fraction < 1.5:
            nice_fraction = 1
        elif fraction < 3:
            nice_fraction = 2
        elif fraction < 7:
            nice_fraction = 5
        else:
            nice_fraction = 10
    else:
        if fraction <= 1:
            nice_fraction = 1
        elif fraction <= 2:
            nice_fraction = 2
        elif fraction <= 5:
            nice_fraction = 5
        else:
            nice_fraction = 10

    return nice_fraction * 10**exponent

def calculate_axis_limits_and_ticks(data_min, data_max, num_ticks=5, exact=False):
    range_span = get_nice_number(data_max - data_min, round=True)
    tick_spacing = get_nice_number(range_span / (num_ticks - 1), round=True)
    if exact:
        nice_min = data_min
        nice_max = data_max
    else:
        nice_min = np.floor(data_min / tick_spacing) * tick_spacing
        nice_max = np.ceil(data_max / tick_spacing) * tick_spacing
    offset = 0.025 * tick_spacing

    return nice_min, nice_max, tick_spacing, offset

def calc_object_viewer_vectors(dataset, current_ind, current_view, viewer_dist):
    # Calculate the viewer-origin vector
    r_v = calc_viewer_position(
        current_view["elev"], current_view["azim"], viewer_dist.to(u.m)
    )

    # Calculate the object-origin vector
    r_o = (
        np.array(
            [
                dataset["x"][current_ind],
                dataset["y"][current_ind],
                dataset["z"][current_ind],
            ]
        )*u.m
    )

    # Calculate the object-viewer vector
    r_ov = r_o - r_v
    return r_v, r_o, r_ov

def calc_object_viewer_angle(r_o, r_ov):
    val = -np.dot(r_ov, r_o) / (np.linalg.norm(r_ov) * np.linalg.norm(r_o))
    object_viewer_angle = np.arccos(
        val
    )
    return object_viewer_angle

def calc_object_viewer_orth_dist(r_v, r_o):
    # Unit vector in the direction of r_vs
    r_v_hat = r_v / np.linalg.norm(r_v)

    # Scalar projection of r_o onto r_v
    scalar_proj = np.dot(r_o, r_v_hat)
    o_v_orth_dist = np.linalg.norm(r_v) - scalar_proj
    return o_v_orth_dist

def create_title(title_dict, animation_key):
    title = ''
    for key, val in title_dict.items():
        if key == animation_key:
            continue
        else:
            title += f"{key}={val}, "
    return title[:-2]
