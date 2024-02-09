import numpy as np


class PlotRenderState:
    def __init__(self, required_keys, base_levels, plot, parallel_map=None):
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
        # Sort the levels of the required keys based on the base_levels
        for level, level_keys in base_levels.items():
            for level_key in level_keys:
                if level_key in self.required_keys:
                    self.levels[level] = level_key
                    self.key_levels[level_key] = level
        # Get the initial values of the coordinates
        # self.next_frame_values = {}
        initial_coords = plot.data.coords
        # for key in self.required_keys:
        #     first_val = initial_coords[key].values[0]
        #     self.next_frame_values[key] = first_val
        # self.last_frame_values = self.next_frame_values.copy()
        self.finished = False
        key_values = {key: initial_coords[key].values for key in self.required_keys}
        self.sel_calls = list(self.generate_all_frames(key_values))
        self.next_sel = self.sel_calls[0]

        if parallel_map is not None:
            # Parallel state is used to keep track of the state of other plots
            # that are running in parallel with this one.
            # parallel_map looks like
            # {"time": "wavelength"}
            # which means the 1st time value is drawn on the same
            # frame as the first "wavelength" value.

            # self.parallel_state = parallel_state
            # parallel_state.parallel_state = self

            # parallel_map is used to map the required keys of this state
            # to the required keys of the parallel state.
            self.parallel_map = {}

            # parallel_values is used to keep track of the values of the current
            # state that need to be compared to the parallel state.
            self.parallel_inds = {}
            for key1, key2 in parallel_map.items():
                if key1 in self.required_keys:
                    self.parallel_inds[key1] = 0
                    self.parallel_map[key1] = key2
                elif key2 in self.required_keys:
                    self.parallel_inds[key2] = 0
                    self.parallel_map[key2] = key1
                else:
                    raise ValueError(
                        f"Neither parallel_map {key1}:{key2} in {self.required_keys}"
                    )

    def iterate_sel(self):
        """
        Updates the next frame values based on the provided new values.

        Args:
            coords (xarray Coordinates):
                A dictionary of the new values for the next frame.
        """
        if self.next_sel == self.sel_calls[-1]:
            self.finished = True
        else:
            self.next_sel = self.sel_calls[self.sel_calls.index(self.next_sel) + 1]
        # keys_updated = False
        # coords = self.plot.data.coords
        # current_key_inds = {}
        # at_last_inds = {}
        # rev_sorted_keys = sorted(
        #     self.required_keys, key=lambda x: self.key_levels[x], reverse=True
        # )
        # for key in rev_sorted_keys:
        #     current_value = self.next_frame_values[key]
        #     key_values = coords[key].values

        #     current_index = np.argwhere(key_values == current_value)[0][0]
        #     current_key_inds[key] = current_index

        #     if current_index == len(key_values) - 1:
        #         at_last_inds[key] = True
        #     else:
        #         at_last_inds[key] = False
        #     # current_ind = np.argwhere(key_values == current_value)
        #     # if current_value == key_values[-1]:
        #     #     # If the current value is the last value, then we're done
        #     #     # with this key on this loop
        #     #     pass
        #     # else:
        #     #     # Otherwise, we need to update the value
        #     #     self.next_frame_values[key] = key_values[current_index + 1]
        #     #     keys_updated = True
        # breakpoint()
        # if all(at_last_inds.values()):
        #     self.finished = True
        # else:
        #     last_key_reset = False
        #     for i, key in enumerate(rev_sorted_keys):
        #         key_values = coords[key].values
        #         if at_last_inds[key]:
        #             last_key_reset = True
        #             self.next_frame_values[key] = key_values[0]

        # if not keys_updated:
        #     # No keys were updated, so we're done
        #     self.finished = True

    def generate_all_frames(self, values_dict):
        # Sort keys by their level
        sorted_keys = sorted(self.key_levels, key=self.key_levels.get)

        def _generate(current_dict, idx):
            if idx == len(sorted_keys):
                yield current_dict.copy()  # We've built a complete combination
            else:
                key = sorted_keys[idx]
                for value in values_dict[key]:
                    current_dict[key] = value
                    yield from _generate(current_dict, idx + 1)

        return _generate({}, 0)

    def is_parent(self, other):
        """
        Check if the other PlotRenderState is the parent of the self PlotRenderState.
        For example, if
            base_levels = {0: ["time"], 1: ["frame"]}
            self.state = {"time": 5, "frame": 1}
            other.state = {"time": 5}
        then other is a parent of self.

        Args:
            other (RenderState):
                The other RenderState instance to compare to.

        Returns:
            bool:
                True if other is a parent state of self
        """
        if self.finished or other.finished:
            return False

        # Check if any of the required keys match
        intersecting_keys = set(self.required_keys).intersection(other.required_keys)
        if not intersecting_keys:
            # No required keys match, so other is not a parent
            return False

        for _, key in sorted(self.levels.items(), key=lambda x: x[0]):
            if key in self.required_keys:
                otherval = other.next_sel.get(key)
                selfval = self.next_sel.get(key)
                if selfval != otherval:
                    # If any key doesn't match, then other is not a parent
                    return False
        return True

    def __lt__(self, other):
        """
        Less than comparison for sorting which RenderState instance will be drawn
        next and is therefore "less".

        Args:
            other (RenderState):
                The other RenderState instance to compare to.

        Returns:
            bool:
                True if self is considered less than other, based on the
                defined key levels and values.
        """
        if self.finished and not other.finished:
            # 'Finished' state is considered greater than any other state
            return False
        if not self.finished and other.finished:
            return True

        for _, key in sorted(self.levels.items(), key=lambda x: x[0]):
            # Keys are checked in levels, keep looking until we find a difference
            # between the two states
            if key in self.required_keys and key in other.required_keys:
                selfval = self.next_sel.get(key)
                otherval = other.next_sel.get(key)
                if selfval != otherval:
                    return selfval < otherval
        return False

    def __eq__(self, other):
        """
        Equal to comparison for checking if two PlotRenderState instances are equal.

        Returns:
            bool:
                True if all required keys and their next frame values are equal
                for both instances.
        """
        if self.finished and other.finished:
            return True
        if self.finished or other.finished:
            return False

        relevant_keys = set(self.required_keys).union(other.required_keys)
        return all(
            self.next_sel.get(key) == other.next_sel.get(key) for key in relevant_keys
        )

    def __repr__(self):
        state_repr = ", ".join(f"{key}={value}" for key, value in self.next_sel.items())
        return f"{self.__class__.__name__}({state_repr}, finished={self.finished})"


class AnimationRenderState:
    def __init__(self, levels):
        """
        The AnimationRenderState class represents the state of the animation
        in full, tracking all values used to render the animation. It allows
        for comparisons between PlotRenderStates that have to do with the state
        of the animation as a whole.

        Args:
            required_keys (list): List of keys this state requires.
            levels (dict): A dictionary mapping keys to their priority level.
            initial_values (dict): Initial values for each key.
        """
        self.levels = levels
        self.animated_values = {}
        self.current_level_inds = {}
        # self.used_values = {key: [] for key in self.levels.keys()}
        for level, level_keys in self.levels.items():
            for level_key in level_keys:
                self.animated_values[level_key] = []
                self.current_level_inds[level_key] = 0

        self.finished = False
