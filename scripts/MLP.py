#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
os.environ["TOKENIZERS_PARALLELISM"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import math
import re
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler, StandardScaler

from transformers import AutoTokenizer
from datasets import Dataset

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader

import lightning as L
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from lightning.pytorch.loggers import TensorBoardLogger, CSVLogger
import torchmetrics
from torchmetrics import Metric


# In[ ]:


print("reading parquet file into datasets...")
df = pd.read_parquet("/home/dzhou/nanopore/parquets/0.5_result_df.parquet")
df.rename(columns={'TRUE_VARIANT': 'labels'}, inplace=True)


###TESTING ONLY
#df = df.sample(frac=0.01)

# Shortening context seq from 201 bp to 101 bp

df['SEQ_CONTEXT'] = df['SEQ_CONTEXT'].apply(lambda x: x[50:151])

# Recalculating GC from shorter seq context

df['GC'] = df['SEQ_CONTEXT'].apply(lambda x : (x.count('G') + x.count('C')) / len(x))

# Adding HOMO field (length of longest homopolymer in the sequence)

def homopolymer_content(seq):
    pattern = r'(A{1,}|C{1,}|T{1,}|G{1,})'
    homopolymers = re.findall(pattern, seq)
    return len(max(homopolymers, key=len))

df['HOMO'] = df['SEQ_CONTEXT'].apply(homopolymer_content)

# Removing DP and SB outliers
df = df[(df['DP'] < 127) & (df['SB'] < 12) & (df['HOMO'] < 12)]
print(df.shape)

# Scaling QUAL, DP, HOMO, and SB columns due to large values

unscaled_numerics = np.transpose(np.stack([df['DP'].values, df['SB'].values, df['HOMO'].values]))
scaler = StandardScaler()
scaled_numerics = scaler.fit_transform(unscaled_numerics)

df['DP'] = scaled_numerics[:,0]
df['SB'] = scaled_numerics[:,1]
df['HOMO'] = scaled_numerics[:,2]

train_val_df, test_df = train_test_split(df, test_size=0.2, stratify=df["labels"])
train_df, val_df = train_test_split(train_val_df, test_size=0.2, stratify=train_val_df["labels"])

train_ds = Dataset.from_pandas(df=train_df, split='train', preserve_index=False)
val_ds = Dataset.from_pandas(df=val_df, split='val', preserve_index=False)
test_ds = Dataset.from_pandas(df=test_df, split='test', preserve_index=False)


# In[ ]:


## HYPERPARAMETERS

tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)

# for model
input_size = 12
hidden_size = 64
num_output_nodes = 1   # for binary classification
dropout = 0.0

# for training/optim

learning_rate = 0.00001
batch_size = 10000
num_epochs = 200


# In[ ]:


def collate_fn(data):
    inputs = [sample['input_ids'] for sample in data]
    labels = torch.tensor([float(sample['labels']) for sample in data]).unsqueeze(1)
    #af = torch.tensor([float(sample['AF']) for sample in data]).unsqueeze(1)    
    #qual = torch.tensor([float(sample['QUAL']) for sample in data]).unsqueeze(1)
    dp = torch.tensor([float(sample['DP']) for sample in data]).unsqueeze(1)
    sb = torch.tensor([float(sample['SB']) for sample in data]).unsqueeze(1)
    homo = torch.tensor([float(sample['HOMO']) for sample in data]).unsqueeze(1)
    gc = torch.tensor([float(sample['GC']) for sample in data]).unsqueeze(1)
    ref = torch.stack([sample['REF'] for sample in data])
    alt = torch.stack([sample['ALT'] for sample in data])
    
    padded_inputs = nn.utils.rnn.pad_sequence(inputs, batch_first=True, padding_value=3) # padding token is 3
    mlp_input = torch.cat((ref, alt, dp, sb, gc, homo), dim=1)  # concatenate other features for mlp input

    return {
        'input_ids': padded_inputs,
        'labels': labels,
        'mlp_input': mlp_input
    }


# In[ ]:

def one_hot_encode(vars, state):
    nuc_dict = {
        'A':np.array([1,0,0,0]),
        'T':np.array([0,1,0,0]),
        'C':np.array([0,0,1,0]),
        'G':np.array([0,0,0,1])
    }
    encoded_vars = []
    for var in vars:
        encoded_vars.append(nuc_dict[var])
    if state == 'REF':
        return {'REF' : encoded_vars}
    elif state == 'ALT':
        return {'ALT' : encoded_vars}
    


# In[ ]:


def dataset_to_dataloader(ds, tokenizer, batch_size=256, shuffle=True, num_workers=4):
    ds = ds.map(lambda x: tokenizer(x['SEQ_CONTEXT'], truncation=False, padding='longest'), batched=True)
    ds = ds.map(lambda x: one_hot_encode(x['REF'], 'REF'), batched=True)   # perform one hot encoding on ref and alt
    ds = ds.map(lambda x: one_hot_encode(x['ALT'], 'ALT'), batched=True)
    ds.set_format(type='torch', columns=['REF', 'ALT', 'DP', 'HOMO', 'GC', 'SB', 'input_ids', 'labels'])
    dl = DataLoader(ds, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers, pin_memory=True, collate_fn=collate_fn)
    return dl


# In[ ]:


print("initializing dataloaders from datasets...")
train_dl = dataset_to_dataloader(train_ds, tokenizer, batch_size=batch_size, shuffle=True)
val_dl = dataset_to_dataloader(val_ds, tokenizer, batch_size=batch_size, shuffle=False)
test_dl = dataset_to_dataloader(test_ds, tokenizer, batch_size=batch_size, shuffle=False)


# In[ ]:


class MLP(L.LightningModule):
    def __init__(self, input_size, hidden_size, num_output_nodes, dropout, lr):
        super(MLP,self).__init__()

        self.dropout = nn.Dropout(dropout)
        self.lr = lr

        self.fc1 = nn.Linear(input_size, hidden_size)
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.fc3 = nn.Linear(hidden_size, num_output_nodes)

        self.train_acc = torchmetrics.classification.BinaryAccuracy()
        self.val_acc = torchmetrics.classification.BinaryAccuracy()
        self.test_acc = torchmetrics.classification.BinaryAccuracy()
        self.test_f1 = torchmetrics.classification.BinaryF1Score()
        self.test_recall = torchmetrics.classification.BinaryRecall()
        self.test_precision = torchmetrics.classification.BinaryPrecision()
        
        self.test_confusion = torchmetrics.classification.BinaryConfusionMatrix()
    
    def forward(self, x):
        x = self.fc1(x)
        x = self.dropout(torch.relu(x))
        x = self.fc2(x)
        x = self.dropout(torch.relu(x))
        x = self.fc3(x)
        return x

    def training_step(self, batch, batch_idx):
        logits = self.forward(batch['mlp_input'])
        targets = batch['labels']

        loss = F.binary_cross_entropy_with_logits(logits, targets)
        scores = torch.sigmoid(logits)
        
        self.train_acc(scores, targets)
        
        self.log_dict(
            {'train_loss':loss, 'train_acc':self.train_acc},
            on_step=False, 
            on_epoch=True, 
            prog_bar=True,
            logger=True
        )
        return loss

    def validation_step(self, batch, batch_idx):
        logits = self.forward(batch['mlp_input'])
        targets = batch['labels']

        loss = F.binary_cross_entropy_with_logits(logits, targets)
        scores = torch.sigmoid(logits)
        
        self.val_acc(scores, targets) 
        
        self.log_dict(
            {'val_loss':loss, 'val_acc':self.val_acc},
            on_step=False,
            on_epoch=True,
            prog_bar=True,
            logger=True
        )
        return loss

    def test_step(self, batch, batch_idx):
        logits = self.forward(batch['mlp_input'])
        targets = batch['labels']

        loss = F.binary_cross_entropy_with_logits(logits, targets)
        scores = torch.sigmoid(logits)
        
        self.test_acc(scores, targets)
        self.test_f1(scores, targets)
        self.test_recall(scores, targets)
        self.test_precision(scores, targets)

        self.test_confusion(scores, targets)

        self.log_dict(
            {'test_loss':loss,
             'test_acc':self.test_acc,
             'test_f1':self.test_f1,
             'test_recall':self.test_recall,
             'test_precision':self.test_precision
            },
            on_step=False,
            on_epoch=True,
            prog_bar=True,
            logger=True
        )
        return loss

    def on_test_epoch_end(self):
        print('CONFUSION_MATRIX:', self.test_confusion.compute())

    
    def configure_optimizers(self):
        optimizer = optim.AdamW(self.parameters(), lr=self.lr)
        return {
        "optimizer": optimizer,
        "lr_scheduler": {
            "scheduler": optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.3, patience=30, threshold=0.001),
            "monitor": "val_acc"
            },
        }


# In[ ]:


model = MLP(
    input_size=input_size,
    hidden_size=hidden_size,
    num_output_nodes=num_output_nodes,
    dropout=dropout,
    lr=learning_rate
)

#logger = TensorBoardLogger('tb_logs', name='mlp_logs')
logger = CSVLogger('csv_logs', name='mlp_csv_logs')

checkpoint_callback = ModelCheckpoint(
    dirpath='/home/dzhou/mlp_checkpoints',
    filename='mlp-{epoch:02d}-{val_acc:.2f}',
    save_top_k=3,  # Save the top 3 models based on monitor
    monitor='val_acc',  # Monitor validation acc
    mode='max',  # maximize validation acc
    verbose=True
)


torch.set_float32_matmul_precision('medium')

trainer = L.Trainer(
    devices=1,
    accelerator="auto",
    max_epochs=num_epochs,
    min_epochs=num_epochs,
    precision="16-mixed",
    logger=logger,
    callbacks=[checkpoint_callback] #EarlyStopping(monitor='val_acc', min_delta=0.0, patience=3)
    #accumulate_grad_batches=2
    #limit_val_batches=0.25,
)

trainer.fit(model, train_dl, val_dl)
trainer.test(model, dataloaders=test_dl)

