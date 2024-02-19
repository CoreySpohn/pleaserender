import numpy as np


class PlotRenderState:
    def __init__(self, required_keys, base_levels, plot, key_strategies=None):
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

        # Get all values for the context that match the required keys
        _context = self.extract_plot_context(fig_context)

        # Change the context and mark for redraw
        self.context = _context
        self.redraw = True

    def extract_plot_context(self, fig_context):
        plot_context = {key: fig_context[key] for key in self.required_keys}
        return plot_context

    # def create_full_context(self, fig_context):
    #     base_context = self.base_sel.copy()
    #     plot_context = self.extract_plot_context(fig_context)
    #     base_context.update(plot_context)
    #     return base_context

    # def check_valid_context(self, fig_context):
    #     # _context = self.extract_plot_context(fig_context)
    #     is_valid = self.plot.valid_plot_data(_context)
    #     return is_valid

    def is_new_context(self, fig_context):
        _context = self.extract_plot_context(fig_context)
        new_context = _context != self.context
        return new_context

    def context_sel(self, context=None):
        if context is None:
            context = self.context
        context = self.extract_plot_context(context)

        data = self.plot.data
        sel = self.base_sel.copy()
        for key, val in context.items():
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
        # sel = self.context_sel(context)
        sel = self.create_full_context(context)
        data = self.plot.data.sel(**sel)
        # for key in self.sum_keys:
        #     if len(data[key].shape) > 0:
        #         data = data.sum(key)

        # data = self.plot.access_data(data)
        return data

    # def is_parent(self, other):
    #     """
    #     Check if the other PlotRenderState is the parent of the self PlotRenderState.
    #     For example, if
    #         base_levels = {0: ["time"], 1: ["frame"]}
    #         self.state = {"time": 5, "frame": 1}
    #         other.state = {"time": 5}
    #     then other is a parent of self.

    #     Args:
    #         other (RenderState):
    #             The other RenderState instance to compare to.

    #     Returns:
    #         bool:
    #             True if other is a parent state of self
    #     """
    #     if self.finished or other.finished:
    #         return False

    #     # Check if any of the required keys match
    #     intersecting_keys = set(self.required_keys).intersection(other.required_keys)
    #     if not intersecting_keys:
    #         # No required keys match, so other is not a parent
    #         return False

    #     for _, key in sorted(self.levels.items(), key=lambda x: x[0]):
    #         if key in self.required_keys:
    #             otherval = other.next_sel.get(key)
    #             selfval = self.next_sel.get(key)
    #             if selfval != otherval:
    #                 # If any key doesn't match, then other is not a parent
    #                 return False
    #     return True

    # def __lt__(self, other):
    #     """
    #     Less than comparison for sorting which RenderState instance will be drawn
    #     next and is therefore "less".

    #     Args:
    #         other (RenderState):
    #             The other RenderState instance to compare to.

    #     Returns:
    #         bool:
    #             True if self is considered less than other, based on the
    #             defined key levels and values.
    #     """
    #     if self.finished and not other.finished:
    #         # 'Finished' state is considered greater than any other state
    #         return False
    #     if not self.finished and other.finished:
    #         return True

    #     for _, key in sorted(self.levels.items(), key=lambda x: x[0]):
    #         # Keys are checked in levels, keep looking until we find a difference
    #         # between the two states
    #         if key in self.required_keys and key in other.required_keys:
    #             selfval = self.next_sel.get(key)
    #             otherval = other.next_sel.get(key)
    #             if selfval != otherval:
    #                 return selfval < otherval
    #     return False

    # def __eq__(self, other):
    #     """
    #     Equal to comparison for checking if two PlotRenderState instances are equal.

    #     Returns:
    #         bool:
    #             True if all required keys and their next frame values are equal
    #             for both instances.
    #     """
    #     if self.finished and other.finished:
    #         return True
    #     if self.finished or other.finished:
    #         return False

    #     relevant_keys = set(self.required_keys).union(other.required_keys)
    #     return all(
    #         self.next_sel.get(key) == other.next_sel.get(key) for key in relevant_keys
    #     )

    def __repr__(self):
        state_repr = ", ".join(f"{key}={value}" for key, value in self.context.items())
        return f"{self.__class__.__name__}({state_repr}, finished={self.finished})"
