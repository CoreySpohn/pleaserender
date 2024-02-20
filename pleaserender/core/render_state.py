import numpy as np


class PlotRenderState:
    def __init__(self, required_keys, base_levels, plot, key_strategies=None, draw_with=None):
        """
        The PlotRenderState class represents the state of a Plot in its render
        process and handles comparisons between different PlotRenderStates.

        Args:
            required_keys (list): List of keys this state requires.
            levels (dict): A dictionary mapping keys to their priority level.
            initial_values (dict): Initial values for each key.
        """
        self.required_keys = required_keys
        self.plot = plot
        self.levels = {}
        self.key_levels = {}
        self.base_sel = plot.base_sel
        self.draw_with = draw_with

        # Sort the levels of the required keys based on the base_levels
        for level, level_keys in base_levels.items():
            for level_key in level_keys:
                if level_key in self.required_keys:
                    self.levels[level] = level_key
                    self.key_levels[level_key] = level

        # Get the initial values of the coordinates
        initial_coords = plot.data.coords
        self.finished = False
        self.key_values = {
            key: initial_coords[key].values for key in self.required_keys
        }

        self.context = None

        # A "key strategy" is a method for handling the key values in the context
        # "Cumulative" means that the key value is used as a slice up to the value
        # "Value" means that the key value is used as a direct value
        # "Trailing" means that the key value is the last value in a 10 value slice
        default_key_strategies = {key: "Value" for key in self.required_keys}
        self.key_strategies = default_key_strategies
        if key_strategies is not None:
            self.key_strategies.update(key_strategies)

    def process_context(self, fig_context):

        # Check if the context is a valid combination for the plot
        valid_context = self.plot.check_valid_context(fig_context)
        if not valid_context:
            # Cannot use this context as the key values do not exist for the plot
            self.redraw = False
            return

        # Check if the context is a new combination
        new = self.is_new_context(fig_context)
        if not new:
            self.redraw = False
            return
        
        if self.draw_with is not None:
            # If this plot is drawn with another plot, then we need to check if
            # the other plot has been redrawn
            if not self.draw_with.state.redraw:
                self.redraw = False
                return

        # Get all values for the context that match the required keys
        _context = self.extract_plot_context(fig_context)

        # Change the context and mark for redraw
        self.context = _context
        self.redraw = True

    def extract_plot_context(self, fig_context):
        plot_context = {key: fig_context[key] for key in self.required_keys}
        return plot_context

    def is_new_context(self, fig_context):
        _context = self.extract_plot_context(fig_context)
        new_context = _context != self.context
        return new_context

    def context_sel(self, context=None, extra_sel=None, extra_strategies=None):
        if context is None:
            context = self.context
        context = self.extract_plot_context(context)

        data = self.plot.data
        sel = self.base_sel.copy()
        if extra_sel is not None:
            # Update the base selection with the extra selection if it exists
            sel.update(extra_sel)
        for key, val in context.items():
            if extra_strategies is not None:
                # If the key has a specific strategy, use that otherwise, use
                # the default strategy
                key_strat = extra_strategies.get(key, self.key_strategies.get(key))
            else:
                key_strat = self.key_strategies.get(key)
            match key_strat:
                case "Cumulative":
                    sel[key] = slice(None, val)
                case "Value":
                    sel[key] = val
                case "Trailing":
                    index = np.where(data[key].values == val)[0][0]
                    sel[key] = slice(max(0, index - 9), index + 1)
                case "All":
                    sel[key] = slice(None, None)
                case _:
                    # Default to "Value"
                    sel[key] = val
        return sel

    def context_data(self, context=None):
        sel = self.create_full_context(context)
        data = self.plot.data.sel(**sel)
        return data

    def __repr__(self):
        state_repr = ", ".join(f"{key}={value}" for key, value in self.context.items())
        return f"{self.__class__.__name__}({state_repr}, finished={self.finished})"
