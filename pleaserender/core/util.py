import numpy as np


def filter_data(dataset, animation_value, animation_key, animation_style):
    match animation_style:
        case "Cumulative":
            return dataset.where(
                    dataset[animation_key] <= animation_value, drop=True)
        case "Single point":
            return dataset.where(
                    dataset[animation_key] == animation_value, drop=True)
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
