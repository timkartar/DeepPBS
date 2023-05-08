import numpy as np
import sys, os
import random

#has_dna = [pdb.strip() for pdb in (open("../extras/has_dna.txt","r").readlines())]

home = "/home/raktim"
data = np.load(home + "/deeppbs/data/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy", allow_pickle = True)

data = [item for item in data if item != []]
random.shuffle(data)
print(len(data))
#correctly_processed = [l.strip() for l in open("/project/rohs_102/share/pwm_data/utils/dataset.txt",'r').readlines()]
correctly_processed = [filename for filename in os.listdir(home + "/deeppbs/data/assembly/") if filename.endswith("npz")]
samples_per_cluster = 5
nclusters_train = 150
score_threshold = 0
align_len_thresold = 4

#prefix = sys.argv[1]

train = []
valid = []

#count = 0
#for item in data:
#    if(item != []):
#        count += 1
#print(count)
#
#sys.exit()
train_clusters = 0
for item in data:
    item_filtered = []
    scores = []
    align_lens = []
    contacts = []
    if(item != []):
        for xy in item:
            xy[1] = np.random.choice(xy[1])
            xy = "_".join(xy) + ".npz"
            if  xy in correctly_processed:
                npz = np.load(home + "/deeppbs/data/assembly/" + xy)
                if npz['contacts'][0] < 50:
                    continue
                item_filtered.append(xy)
                scores.append(npz['aln_score'][0])
                align_lens.append(npz['pwm_mask'][0].sum())
                contacts.append(npz['contacts'][0])
        
        #arg = np.argsort(scores)[::-1]
        #arg = list(range(len(scores)))

        ### filter by contacts with aligned region
        arg = np.argsort(contacts)[::-1]  

        #np.random.shuffle(arg)
        
        item_filtered = np.array(item_filtered)[arg]
        scores = np.array(scores)[arg]
        align_lens = np.array(align_lens)[arg]
        count = 0

        if(train_clusters < nclusters_train):
            for i in range(len(item_filtered)):

                #if(xy[0] not in has_dna):
                toappend = item_filtered[i]
                if scores[i] > score_threshold and align_lens[i] > align_len_thresold:
                    train.append(item_filtered[i])
                    count += 1
                if(count >= samples_per_cluster):
                    break
            if(count > 0):
                train_clusters += 1
        else:
            for i in range(len(item_filtered)):
                #if(xy[0] not in has_dna):
                toappend = item_filtered[i]
                if scores[i] > score_threshold and align_lens[i] > align_len_thresold:
                    valid.append(item_filtered[i])
                    count += 1
                if(count >= samples_per_cluster):
                    break
    


#open("./train0" + str(samples_per_cluster)  + ".txt","w").write("\n".join(train))
#open("./valid0"+ str(samples_per_cluster) + ".txt","w").write("\n".join(valid))
print(len(train), len(valid))

all_data = train + valid


nfolds = 5

foldsize = len(all_data)//nfolds

'''
for i in range(nfolds):
    valid = all_data[i*foldsize:(i+1)*foldsize]
    open("./folds/valid" + str(i)  + ".txt","w").write("\n".join(valid))
    train = [item for item in all_data if item not in valid]
    open("./folds/train" + str(i)  + ".txt","w").write("\n".join(train))
'''
