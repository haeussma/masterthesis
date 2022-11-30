from typing import List

import numpy as np
import pandas as pd
import re
from datetime import datetime

def blockshaped(arr: np.ndarray, nrows: int, ncols: int):
        """
        Return an array of shape (n, nrows, ncols) where
        n * nrows * ncols = arr.size

        If arr is a 2D array, the returned array should look like n subblocks with
        each subblock preserving the "physical" layout of arr.
        """
        h, w = arr.shape
        assert h % nrows == 0, f"{h} rows is not evenly divisible by {nrows}"
        assert w % ncols == 0, f"{w} cols is not evenly divisible by {ncols}"
        return (arr.reshape(h//nrows, nrows, -1, ncols)
                .swapaxes(1,2)
                .reshape(-1, nrows, ncols))


def read_photometer(path) -> pd.DataFrame:
        return pd.read_csv(
                path, 
                sep='delimiter', 
                encoding='utf-16', 
                engine='python')


def string_to_float(string: str) -> float:
        number = re.sub(r'[^0-9.]', '', string)
        if len(number) == 0:
            return float('nan')
        else:
            return float(number)

def tab_split(string: str) -> list:
        return string.split("\t")


def to_seconds(string: str) -> List[int]:
        time = datetime.strptime(string, "%H:%M:%S").time()
        return time.hour*3600 + time.minute*60 + time.second