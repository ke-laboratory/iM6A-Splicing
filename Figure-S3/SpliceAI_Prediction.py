from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import pandas as pd

## Load models
context = 10000
paths = ('/home/morgenlefay/iM6A/SpliceAI/models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(y) for y in paths]
models

## Load data
Fasta = pd.read_csv("/home/morgenlefay/iM6A/Chapter3/Figure2/Splicing/mm10_Fasta1.csv")
Sequence = Fasta["Seq"].tolist()

## Prediction for WT intron
for i in range(len(Fasta)):
    print(i)
    input_sequence = Sequence[i]

    x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
    y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

    donor_prob = y[0, :, 2]
    acceptor_prob = y[0, :, 1]
    donor_prob = donor_prob.tolist()
    acceptor_prob = acceptor_prob.tolist()

    df = pd.DataFrame(np.random.randn((len(input_sequence)), 3))
    df.columns = ["name", "chrom", "strand"]
    df["name"] = Fasta.loc[i,"name"]
    df["chrom"] = Fasta.loc[i,"chrom"]
    df["strand"] = Fasta.loc[i,"strand"]
    df = df.reset_index(drop = False)
    df.rename(columns={"index":"Position"}, inplace=True)

    Probability = pd.DataFrame({'donor_prob':donor_prob,'acceptor_prob':acceptor_prob})
    df = pd.concat([df, Probability], axis=1)
    df = df[(df["donor_prob"]>=0.1) | (df["acceptor_prob"]>=0.1)]
    df = df.reset_index(drop = True)
    df.to_csv("/home/morgenlefay/iM6A/Chapter3/Figure2/Splicing/Output1/{}.bed".format(Fasta.loc[i,"name"]), sep="\t", index=False)    
