class RenderState:
    def __init__(self, required_keys, order, initial_values):
        """
        The RenderState class represents the state of the render process and
        handles comparisons between different states.

        Args:
            required_keys (list): List of keys this state requires.
            order (dict): A dictionary mapping keys to their priority order.
            initial_values (dict): Initial values for each key.
        """
        self.required_keys = required_keys
        self.order = order
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

        for key in sorted(self.order.keys(), key=lambda x: self.order[x]):
            # Keys are checked in order, keep looking until we find a difference
            # between the two states
            if key in self.required_keys and key in other.required_keys:
                selfval = self.next_frame_values.get(key)
                otherval = other.next_frame_values.get(key)
                # return value1 < value2
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

        return all(
            self.next_frame_values.get(key) == other.next_frame_values.get(key)
            for key in self.required_keys
            if key in other.required_keys
        )
