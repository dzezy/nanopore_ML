#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import re
import os
os.environ["TOKENIZERS_PARALLELISM"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import math
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from datasets import Dataset
from transformers import AutoTokenizer, PreTrainedTokenizerFast

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
#df = pd.read_parquet("/home/dzhou/nanopore/parquets/sampled_parsed_result_df.parquet")
df.rename(columns={'TRUE_VARIANT': 'labels'}, inplace=True)


###TESTING ONLY
#df = df.sample(frac=0.001)

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

# Removing DP, HOMO and SB outliers (keeping only 99th percentile of data for each)
df = df[(df['DP'] < 127) & (df['SB'] < 12) & (df['HOMO'] < 12)]
# leaves 10.15 M compared to 10.5M original dataset


# Scaling QUAL, DP, HOMO, and SB columns due to large values

unscaled_numerics = np.transpose(np.stack([df['DP'].values, df['SB'].values, df['HOMO'].values]))
scaler = StandardScaler()
scaled_numerics = scaler.fit_transform(unscaled_numerics)

df['DP'] = scaled_numerics[:,0]
df['SB'] = scaled_numerics[:,1]
df['HOMO'] = scaled_numerics[:,2]


# Creating data splits
train_val_df, test_df = train_test_split(df, test_size=0.2, stratify=df["labels"])
train_df, val_df = train_test_split(train_val_df, test_size=0.2, stratify=train_val_df["labels"])


# In[ ]:


train_ds = Dataset.from_pandas(df=train_df, split='train', preserve_index=False)
val_ds = Dataset.from_pandas(df=val_df, split='val', preserve_index=False)
test_ds = Dataset.from_pandas(df=test_df, split='test', preserve_index=False)


# In[ ]:


tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
'''
tokenizer = PreTrainedTokenizerFast(
    tokenizer_file="/home/dzhou/0.5_dataset_BPEtokenizer.json",
    pad_token = '[PAD]',
    sep_token = '[SEP]',
    unk_token = '[UNK]',
    cls_token = '[CLS]',
    mask_token = '[MASK]'
)
'''

## HYPERPARAMETERS

# for model
src_vocab_size = tokenizer.vocab_size
d_model = 512
num_heads = 8
num_layers = 2
d_ff = 2048
max_seq_len = 35
dropout = 0.4
mean_pool = True

# for training/optim

learning_rate = 0.00003
batch_size = 2048
num_epochs = 100
num_output_nodes = 1   # for binary classification


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


# In[9]:


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

def dataset_to_dataloader(ds, tokenizer, batch_size=256, shuffle=True, num_workers=4):
    ds = ds.map(lambda x: tokenizer(x['SEQ_CONTEXT'], truncation=False, padding='longest'), batched=True)
    ds = ds.map(lambda x: one_hot_encode(x['REF'], 'REF'), batched=True)   # perform one hot encoding on ref and alt
    ds = ds.map(lambda x: one_hot_encode(x['ALT'], 'ALT'), batched=True)
    ds.set_format(type='torch', columns=['REF', 'ALT', 'DP', 'HOMO', 'GC', 'SB', 'input_ids', 'labels'])
    dl = DataLoader(ds, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers, pin_memory=True, collate_fn=collate_fn)
    return dl


print("initializing dataloaders from datasets...")
train_dl = dataset_to_dataloader(train_ds, tokenizer, batch_size=batch_size, shuffle=True)
val_dl = dataset_to_dataloader(val_ds, tokenizer, batch_size=batch_size, shuffle=False)
test_dl = dataset_to_dataloader(test_ds, tokenizer, batch_size=batch_size, shuffle=False)


# In[ ]:


class MultiHeadAttention(nn.Module):
    def __init__(self, d_model, num_heads):
        super(MultiHeadAttention, self).__init__()
        
        # assert model dimension is divisible by number of heads
        assert d_model % num_heads == 0, "d_model must be divisible by num_heads"

        self.d_model = d_model  # model dimension, aka embedding size
        self.num_heads = num_heads   # number of attn heads
        self.d_k = d_model // num_heads    # dim of each head's key, query, and value

        self.W_q = nn.Linear(d_model, d_model) 
        self.W_k = nn.Linear(d_model, d_model) 
        self.W_v = nn.Linear(d_model, d_model) 
        self.W_o = nn.Linear(d_model, d_model)

    def split_heads(self, x):
        # Reshape the input to have num_heads for multi-head attention
        # (batch, seq_len, d_model) --> (batch, num_heads, seq_len, d_k)
        batch_size, seq_length, d_model = x.size()
            
        return x.view(batch_size, seq_length, self.num_heads, self.d_k).transpose(1, 2)

    def combine_heads(self, x):
        # Combine the multiple heads back to original shape
        batch_size, _, seq_length, d_k = x.size()
        return x.transpose(1, 2).contiguous().view(batch_size, seq_length, self.d_model)
    
    
    def forward(self, Q, K, V, mask=None):
        # apply linear transformations and split heads
        Q = self.split_heads(self.W_q(Q))
        K = self.split_heads(self.W_k(K))
        V = self.split_heads(self.W_v(V))

        #with torch.backends.cuda.sdp_kernel(enable_flash=True, enable_math=False, enable_mem_efficient=False):
            # scaled dot-product attention for query key value
        attn_output = F.scaled_dot_product_attention(Q, K, V, attn_mask=mask, dropout_p=0.0)

        # combine heads and apply output transformation
        return self.W_o(self.combine_heads(attn_output))


# In[ ]:


class FeedForward(nn.Module):
    def __init__(self, d_model, d_ff):
        super(FeedForward, self).__init__()

        # initialize feedforward layers and activation 
        self.fc1 = nn.Linear(d_model, d_ff)
        self.fc2 = nn.Linear(d_ff, d_model)
        self.gelu = nn.GELU()  # using GELU instead of ReLU

    def forward(self, x):
        return self.fc2(self.gelu(self.fc1(x)))


# In[ ]:


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_seq_len):
        super(PositionalEncoding, self).__init__()

        pe = torch.zeros(max_seq_len, d_model)
        position = torch.arange(0, max_seq_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * -(math.log(10000.0) / d_model))

        # using sin and cosine formulas for positional encoding
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        # save positional encoding as static buffer (not trainable param) of the model
        self.register_buffer('pe', pe.unsqueeze(0)) 

    def forward(self, x):
        return x + self.pe[:, :x.size(1)]
        


# In[ ]:


class EncoderLayer(nn.Module):
    def __init__(self, d_model, num_heads, d_ff, dropout):
        super(EncoderLayer, self).__init__()
        self.self_attn = MultiHeadAttention(d_model, num_heads)
        self.feed_forward = FeedForward(d_model, d_ff)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, mask):
        
        attn_output = self.self_attn(x, x, x, mask)
        x = self.norm1(x + self.dropout(attn_output)) # residual/skip connections 
        ff_output = self.feed_forward(x)
        x = self.norm2(x + self.dropout(ff_output))
        return x


# In[ ]:


class Transformer(L.LightningModule):
    def __init__(self, src_vocab_size, d_model, num_heads, num_layers, d_ff, max_seq_len, num_output_nodes, dropout=0.0, mean_pool=True, lr=0.00001):
        super(Transformer, self).__init__()

        self.encoder_embedding = nn.Embedding(src_vocab_size, d_model)
        self.positional_encoding = PositionalEncoding(d_model, max_seq_len)
        self.encoder_layers = nn.ModuleList([EncoderLayer(d_model, num_heads, d_ff, dropout) for _ in range(num_layers)])

        self.fc = nn.Linear(d_model, num_output_nodes)
        self.dropout = nn.Dropout(dropout)
        self.mean_pool = mean_pool
        self.lr = lr
        
        self.train_acc = torchmetrics.classification.BinaryAccuracy()
        self.val_acc = torchmetrics.classification.BinaryAccuracy()
        self.test_acc = torchmetrics.classification.BinaryAccuracy()
        self.test_f1 = torchmetrics.classification.BinaryF1Score()
        self.test_recall = torchmetrics.classification.BinaryRecall()
        self.test_precision = torchmetrics.classification.BinaryPrecision()
        self.test_confusion = torchmetrics.classification.BinaryConfusionMatrix()
    
    def forward(self, src):
        src_mask = (src != 3).unsqueeze(1).unsqueeze(2)  # mask out padding tokens (where input_id = 3)
        src_embedded = self.dropout(self.positional_encoding(self.encoder_embedding(src)))
        
        enc_output = src_embedded
        for enc_layer in self.encoder_layers:
            enc_output = enc_layer(enc_output, src_mask)

        # perform mean or max pooling over seq_len dimension
        pooled_output = enc_output.mean(dim=1) if self.mean_pool else enc_output.max(dim=1)[0]
        
        final_output = self.fc(pooled_output)
        return final_output

    def training_step(self, batch, batch_idx):
        logits = self.forward(batch['input_ids'])
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
        logits = self.forward(batch['input_ids'])
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
        logits = self.forward(batch['input_ids'])
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
            }
        }


# In[ ]:


model = Transformer(
    src_vocab_size, 
    d_model, 
    num_heads, 
    num_layers, 
    d_ff, 
    max_seq_len, 
    num_output_nodes, 
    dropout, 
    mean_pool, 
    learning_rate
)


#logger = TensorBoardLogger('tb_logs', name='transformer_logs')
logger = CSVLogger('csv_logs', name='short_transformer_csv_logs')


checkpoint_callback = ModelCheckpoint(
    dirpath='/home/dzhou/short_transformer_checkpoints',
    filename='short_transformer-{epoch:02d}-{val_acc:.2f}',
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
    #accumulate_grad_batches=2
    #limit_val_batches=0.25,
    callbacks=[checkpoint_callback] #, EarlyStopping(monitor='val_acc', min_delta=0.0, patience=30)]
)

trainer.fit(model, train_dl, val_dl)
trainer.test(model, dataloaders=test_dl)

