# import pandas as pd
import xarray as xr


class PlotObject:
    def __init__(self, input):
        """
        Thin wrapper around a plot object, such as a planetary system,
        so that dataframes can be used if the data is not generated.
        """
        if isinstance(input, xr.Dataset):
            self.data = input
            self.is_xr = True
        else:
            self.obj = input
            self.data = xr.Dataset()
            self.is_xr = False

    # def filter_data(self, animation_value, animation_key, animation_style):
    #     match animation_style:
    #         case "Cumulative":
    #             return self.data[self.data[animation_key] <= animation_value]
    #         case "Single point":
    #             return self.data[self.data[animation_key] == animation_value]
    #         case "Trailing":
    #             # Data of trailing 10 points
    #             current_index = (self.data[animation_key] == animation_value).idxmax()
    #             min_index = self.data.index[0]
    #             if current_index - 10 > min_index:
    #                 return self.data[
    #                     (self.data.index <= current_index)
    #                     & (self.data.index >= current_index - 10)
    #                 ]
    #             else:
    #                 return self.data[self.data.index <= current_index]

    def filter_data(self, animation_value, animation_key, animation_style):
        match animation_style:
            case "Cumulative":
                return self.data.where(
                    self.data[animation_key] <= animation_value, drop=True
                )
            case "Single point":
                return self.data.where(
                    self.data[animation_key] == animation_value, drop=True
                )
            case "Trailing":
                # Data of trailing 10 points
                # Finding the index of the specified animation value
                current_index = self.data[animation_key].argmax(dim="your_dimension")
                min_index = 0  # Assuming 0 is the starting index in your dimension

                if current_index - 10 > min_index:
                    return self.data.where(
                        (self.data[animation_key].dims[0] <= current_index)
                        & (self.data[animation_key].dims[0] >= current_index - 10),
                        drop=True,
                    )
                else:
                    return self.data.where(
                        self.data[animation_key].dims[0] <= current_index, drop=True
                    )

    # def get_data_rows(self, keys, axis_values, axis_name="t"):
    #     """
    #     Returns a DataFrame with the data required to render the plot.
    #     Args:
    #         keys (list): List of keys to be returned, such as ['x', 'y']
    #     """
    #     if self.is_df:
    #         return self.please_render_data
