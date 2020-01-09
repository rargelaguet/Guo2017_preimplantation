from biofam.run.entry_point import entry_point
import pandas as pd
from time import time

# infile = "/Users/ricard/data/Guo_2017/mofa/data.txt.gz"
# outfile = "/Users/ricard/data/Guo_2017/mofa/test.hdf5"

infile = "/hps/nobackup2/research/stegle/users/ricard/Guo_2017/mofa/data.txt.gz"
outfile = "/hps/nobackup2/research/stegle/users/ricard/Guo_2017/mofa/hdf5/model_test_large.hdf5"


data = pd.read_csv(infile, delimiter="\t", header=0)
lik = ["gaussian"]*data["feature_group"].nunique()

ent = entry_point()
ent.set_data_options(likelihoods=lik, center_features_per_group=True, scale_features=False, scale_views=False)
ent.set_data_df(data)

ent.set_model_options(factors=25, likelihoods=lik, sl_z=True, sl_w=True, ard_z=False, ard_w=True)

ent.set_train_options(iter=5000, convergence_mode="medium", dropR2=None, gpu_mode=True, startELBO=25, elbofreq=3, verbose=False)

ent.build()
ent.run()

ent.save(outfile)