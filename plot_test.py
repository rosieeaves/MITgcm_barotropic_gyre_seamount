import numpy as np
import xmitgcm as mit

data = mit.open_mdsdataset('./run',iters=0, read_grid=False)
print(data.U.dims)