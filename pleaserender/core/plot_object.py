import pandas as pd


class PlotObject:
    def __init__(self, input):
        """
        Thin wrapper around a plot object, such as a planetary system,
        so that dataframes can be used if the data is not generated.
        """
        if isinstance(input, pd.DataFrame):
            self.data = input
            self.is_df = True
        else:
            self.obj = input
            self.data = pd.DataFrame()
            self.is_df = False

    def filter_data(self, animation_value, animation_key, animation_style):
        match animation_style:
            case "Cumulative":
                return self.data[self.data[animation_key] <= animation_value]
            case "Single point":
                return self.data[self.data[animation_key] == animation_value]
            case "Trailing":
                # Data of trailing 10 points
                current_index = (self.data[animation_key] == animation_value).idxmax()
                min_index = self.data.index[0]
                if current_index - 10 > min_index:
                    return self.data[
                        (self.data.index <= current_index)
                        & (self.data.index >= current_index - 10)
                    ]
                else:
                    return self.data[self.data.index <= current_index]


    # def get_data_rows(self, keys, axis_values, axis_name="t"):
    #     """
    #     Returns a DataFrame with the data required to render the plot.
    #     Args:
    #         keys (list): List of keys to be returned, such as ['x', 'y']
    #     """
    #     if self.is_df:
    #         return self.please_render_data
