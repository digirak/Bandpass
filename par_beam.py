from mwa_new_beam import mwa_beam
from multiprocessing import Process as p
freqs=np.arange(100e+6,200e+6,1e+6)
n=100
pool=p(target=mwa_beam,args=(n,freqs))
pool.start()
pool.join()
