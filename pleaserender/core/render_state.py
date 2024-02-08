class RenderState:
    def __init__(self, required_keys, base_order, initial_values):
        """
        The RenderState class represents the state of the render process and
        handles comparisons between different states.

        Args:
            required_keys (list): List of keys this state requires.
            order (dict): A dictionary mapping keys to their priority order.
            initial_values (dict): Initial values for each key.
        """
        self.required_keys = required_keys
        self.order = {}
        # Sort the order of the required keys based on the base_order
        for order_int, order_vals in base_order.items():
            for order_val in order_vals:
                if order_val in self.required_keys:
                    self.order[order_int] = order_val
        self.last_frame_values = initial_values
        self.next_frame_values = initial_values.copy()
        self.finished = False

    def update_next_frame_values(self, coords):
        """
        Updates the next frame values based on the provided new values.

        Args:
            coords (xarray Coordinates):
                A dictionary of the new values for the next frame.
        """
        keys_updated = False
        for key in sorted(self.required_keys, key=lambda x: self.order[x]):
            current_value = self.next_frame_values[key]
            available_values = coords[key]

            if current_value in available_values:
                current_index = available_values.index(current_value)
                if current_index + 1 < len(available_values):
                    self.next_frame_values[key] = available_values[current_index + 1]
                    keys_updated = True
                    break
                else:
                    self.next_frame_values[key] = available_values[0]
            else:
                breakpoint()

        if not keys_updated:
            # No keys were updated, so we're done
            self.finished = True

    def is_parent(self, other):
        """
        Check if the other RenderState is the parent of the self RenderState.
        For example, if
            base_order = {0: ["time"], 1: ["frame"]}
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

        for _, key in sorted(self.order.items(), key=lambda x: x[0]):
            if key in self.required_keys:
                otherval = other.next_frame_values.get(key)
                selfval = self.next_frame_values.get(key)
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
                defined key order and values.
        """
        if self.finished and not other.finished:
            # 'Finished' state is considered greater than any other state
            return False
        if not self.finished and other.finished:
            return True

        for _, key in sorted(self.order.items(), key=lambda x: x[0]):
            # Keys are checked in order, keep looking until we find a difference
            # between the two states
            if key in self.required_keys and key in other.required_keys:
                selfval = self.next_frame_values.get(key)
                otherval = other.next_frame_values.get(key)
                if selfval != otherval:
                    return selfval < otherval
        return False

    def __eq__(self, other):
        """
        Equal to comparison for checking if two RenderState instances are equal.

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
            self.next_frame_values.get(key) == other.next_frame_values.get(key)
            for key in relevant_keys
        )

    def __repr__(self):
        state_repr = ", ".join(
            f"{key}={value}" for key, value in self.next_frame_values.items()
        )
        return f"{self.__class__.__name__}({state_repr}, finished={self.finished})"
