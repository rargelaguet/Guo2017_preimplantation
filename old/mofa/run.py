from mofapy2.run.entry_point import entry_point
import pandas as pd

infile = "/Users/ricard/data/Guo_2017/mofa/data.txt.gz"
outfile = "/Users/ricard/data/Guo_2017/mofa/hdf5/model.hdf5"

data = pd.read_csv(infile, delimiter="\t", header=0)
lik = ["gaussian"]*data["view"].nunique()

ent = entry_point()
ent.set_data_options(likelihoods=lik)
ent.set_data_df(data)

ent.set_model_options(factors=5, likelihoods=lik, spikeslab_weights=True, ard_factors=False, ard_weights=True)

ent.set_train_options(iter=5000, convergence_mode="fast", startELBO=1, elbofreq=1, verbose=False)

ent.build()
ent.run()

ent.save(outfile)