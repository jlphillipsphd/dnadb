import numpy as np

class MockableRng(np.random.Generator):
    def __init__(self):
        super().__init__(np.random.PCG64())

    def choice(self, *args, **kwargs):
        return super().choice(*args, **kwargs)
