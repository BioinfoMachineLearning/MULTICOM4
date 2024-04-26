from transformers import AutoTokenizer, EsmForMaskedLM, pipeline
import argparse
import os, sys
import numpy as np
from torch.utils.data import Dataset
from tqdm import tqdm
import random

class MyDataset(Dataset):
    def __init__(self, sequence):
        self.sequence = sequence
        self.batched_sequence = []
        self._prepare()

    def __len__(self):
        return len(self.batched_sequence)

    def __getitem__(self, idx):
        return self.batched_sequence[idx]

    def _prepare(self):
        for idx in range(len(self.sequence)):
            residue = self.sequence[idx]
            text = self.sequence[0:idx] + "<mask>" + self.sequence[idx+1:]
            self.batched_sequence += [text]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ina3m', type=str, required=True)
    parser.add_argument('--outfile', type=str, required=True)
    parser.add_argument('--threshold', type=float, required=False, default=0.1)
    args = parser.parse_args()

    a3m_contents = []

    for line in open(args.ina3m):
        a3m_contents += [line.rstrip('\n')]

    sequence = a3m_contents[1]

    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t36_3B_UR50D")
    model = EsmForMaskedLM.from_pretrained("facebook/esm2_t36_3B_UR50D")
    pipe = pipeline("fill-mask", model=model, tokenizer=tokenizer, device=0)

    dataset = MyDataset(sequence)

    mask_predictions = []
    for out in tqdm(pipe(dataset, batch_size=3, top_k=2), total=len(dataset)):
        mask_predictions += [out]

    list_0 = []
    list_1 = []
    for idx in range(len(sequence)):
        mask_pred = mask_predictions[idx]
        if mask_pred[0]["score"] > args.threshold:
            list_0.append(mask_pred[0]["token_str"])
        else:
            list_0.append("-")

        if mask_pred[1]["score"] > args.threshold*0.9:
            list_1.append(mask_pred[1]["token_str"])
        else:
            list_1.append("-")

    sampled_sequence_1 = ''.join(list_0)
    sampled_sequence_2 = ''.join(list_1)

    a3m_contents += [f">alg1_0"]
    a3m_contents += [sampled_sequence_1]
        
    a3m_contents += [f">alg1_1"]
    a3m_contents += [sampled_sequence_2]

    open(args.outfile, 'w').write('\n'.join(a3m_contents))


