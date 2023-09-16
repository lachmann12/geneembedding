import random
import openai
import json
import tqdm
import numpy as np
import pandas as pd

def read_gmt(file):
    with open(file, 'r') as f:
        library = {line.split('\t')[0]: line.strip().split('\t')[2:] for line in f}
    return library

def split_lists(input_dict, factor=1, overlap=False):
    dict1 = {}
    dict2 = {}
    for key, value in input_dict.items():
        for i in range(factor):
            random.shuffle(value)
            split_point = len(value) // 2
            dict1[key+str(i)] = value[:split_point]
            dict2[key+str(i)] = value[split_point:]
            dict1[key+str(i)+"o"] = random.sample(value, split_point)
            dict2[key+str(i)+"o"] = random.sample(value, split_point)
    return dict1, dict2

def chunk_list(lst, chunk_size):
    return [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]

def remove_long_lists(dictionary, n):
    return {key: value for key, value in dictionary.items() if not isinstance(value, list) or len(value) <= n}

def embedding_vec(data, model="text-embedding-ada-002", max_chunk_size=50, progress=False):
    with open("secrets/config.json") as file:
        openai.api_key = json.loads(file.read())["chatgpt"]["key"]
    chunk_data = chunk_list(data, max_chunk_size)
    embeddings = []
    for chunk_data in tqdm.tqdm(chunk_data, disable=not progress):
        embedding_result = openai.Embedding.create(input=chunk_data, model=model)
        chunk_embeddings = np.array([x["embedding"] for x in embedding_result["data"]])
        embeddings.append(chunk_embeddings)
    return np.concatenate(embeddings, axis=0)

def embed_library(library, progress=False):
    concatenated_list = [','.join(value) for key, value in library.items()]
    return pd.DataFrame(embedding_vec(concatenated_list, progress=progress), index=list(library.keys()))