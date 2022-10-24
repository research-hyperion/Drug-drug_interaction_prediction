
import networkx as nx 
import matplotlib 
import matplotlib.pyplot as plt 
import numpy as np
import collections
import powerlaw
import operator
import scipy.stats
import random
import math
import csv
import random
import pandas as pd
from tqdm import tqdm


"""https://networkx.org/documentation/stable/tutorial.html#drawing-graphs

https://networkx.org/documentation/stable/reference/algorithms/link_prediction.html
"""

DRUG_FILE = "Classified_DDis-drugs.com-Sheet1.csv"

df = pd.read_csv(DRUG_FILE)

all_interactions = []
severe_moderate_interactions = []
severe_interactions=[]
all_nodes = []
result_list = []
for index, row in df.iterrows():
    drug_1 = row['DB 1']
    drug_1_name = row['Drug 1']
    drug_2 = row['DB2']
    drug_2_name = row['Drug 2']
    
    new_entry = {"drugname1":drug_1_name, "dbid1":drug_1, "drugname2":drug_2_name, "dbid2":drug_2, "score_low":0, "score_medium":0, "score_high":0}
    result_list.append(new_entry)
    interaction_strength = row["Interaction strength"]
    if interaction_strength == "major":
        severe_interactions.append((drug_1,drug_2))
        severe_moderate_interactions.append((drug_1,drug_2))
        all_interactions.append((drug_1,drug_2))
        continue
    elif interaction_strength == "moderate":
        severe_moderate_interactions.append((drug_1,drug_2))
        all_interactions.append((drug_1,drug_2))
        continue
    elif interaction_strength == "minor":
        all_interactions.append((drug_1,drug_2))
        continue
    else:
        all_nodes.append(drug_1)
        all_nodes.append(drug_2)

all_interactions_graph = nx.Graph()
severe_moderate_graph = nx.Graph()
severe_interactions_graph = nx.Graph()

all_interactions_graph.add_edges_from(all_interactions)
severe_moderate_graph.add_edges_from(severe_moderate_interactions)
severe_interactions_graph.add_edges_from(severe_interactions)

all_interactions_graph.add_nodes_from(all_nodes)
severe_moderate_graph.add_nodes_from(all_nodes)
severe_interactions_graph.add_nodes_from(all_nodes)

#All interactions
non_edges = list(nx.non_edges(all_interactions_graph))
preds_all_interactions = nx.adamic_adar_index(all_interactions_graph, non_edges)


# for u, v, p in preds_all_interactions:
#     print(f"({u}, {v}) -> {p:.8f}")

non_edges = list(nx.non_edges(severe_moderate_graph))
preds_severe_moderate_interactions = nx.adamic_adar_index(severe_moderate_graph, non_edges)

non_edges = list(nx.non_edges(severe_interactions_graph))
preds_severe_interactions = nx.adamic_adar_index(severe_interactions_graph, non_edges)

from threading import Thread, Lock
import threading
mutex = Lock()

#{"drugname1":drug_1_name, "dbid1":drug_1, "drugname2":drug_2_name, "dbid2":drug_2, "score_low":0, "score_medium":0, "score_high":0}
def search_in_list(id_in_list,search_list, tag):

    entry = result_list[id_in_list]

    for u,v,p in search_list:
        if entry['dbid1'] == u and entry['dbid2']==v:
            entry[tag] = p
            break
    mutex.acquire()
    result_list[id_in_list]=entry
    mutex.release()

for i in tqdm(range(len(result_list))):
     

    threads = list()
    threads.append(threading.Thread(target=search_in_list, args=(i,preds_all_interactions,"score_low")))
    threads.append(threading.Thread(target=search_in_list, args=(i,preds_severe_moderate_interactions,"score_medium")))
    threads.append(threading.Thread(target=search_in_list, args=(i,preds_severe_interactions,"score_high")))

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

with open('all_intearactions_scores.csv', 'w', newline='') as csvfile:
    fieldnames = ['Drugname1','DBid1', 'Drugname2','DBid2', 'Score_low', 'Score_medium', 'Score_high']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for entry in result_list:
        writer.writerow({'Drugname1': entry['drugname1'], 
                         'DBid1':entry['dbid1'],
                         'Drugname2':entry["drugname2"],
                         'DBid2': entry['dbid2'],
                         'Score_low':entry['score_low'],
                         'Score_medium':entry['score_medium'],
                         'Score_high':entry['score_high']
                         })

