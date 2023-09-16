import pandas as pd
import numpy as np

eranks = pd.DataFrame()

def rel_score(enrichment_ranks, relevance_ranks):
    # normalize ranks to numbers between 0-100