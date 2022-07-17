#!/usr/bin/env python
# coding: utf-8

import math
import random
import csv
import os, time
import re
from multiprocessing import Process, Pool
import pandas as pd
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, MaxPool1D, Flatten, Dropout, Conv1D
from tensorflow.keras.callbacks import EarlyStopping, TensorBoard, ModelCheckpoint
from tensorflow.keras import regularizers
import tensorflow as tf
from sklearn.metrics import confusion_matrix

#parse fasta file into name, seq
def get_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))
def get_1000_window(input_bed, genome_dir, out_dir):
    chrm = []
    start = [] 
    end = []
    with open(input_bed) as in_f:
        reader = csv.reader(in_f, delimiter='\t')
        for row in reader:
            chrm.append(row[0])
            start.append(row[1])
            end.append(row[2])
    chrm = np.array(chrm)
    start = np.array(start, dtype = int)
    end = np.array(end, dtype = int)
    bed = np.stack([chrm, start, end]).transpose()
    print(f"Number of lines before processing: {len(chrm)}")
    #obtain 1000bp window centered on eccDNA midpoint

    smaller_1000 = bed[((bed[:, 2].astype(int) - bed[:, 1].astype(int)) <= 1000), :]
    larger_1000 = bed[((bed[:, 2].astype(int) - bed[:, 1].astype(int)) > 1000), :]
    window_1000 = bed[((bed[:, 2].astype(int) - bed[:, 1].astype(int)) <= 1000), :]
    shift = (1000 - (window_1000[:,2].astype(int) - window_1000[:,1].astype(int)))/2
    window_1000[:, 1] = (window_1000[:, 1].astype(int) - np.floor(shift)).astype(int)
    window_1000[:, 2] = (window_1000[:, 2].astype(int) + np.ceil(shift)).astype(int)
    print("Number of records with 1000 bp window:", len(window_1000))

    hg38 = {}
    with open(genome_dir) as in_g:
        reader = csv.reader(in_g, delimiter='\t')
        for row in reader:
            hg38[row[0]] = row[1]

    #adjust coordinate if < 0 or > chromosome size
    #if coordinate < 0 then change to 0
    lower_than_0 = np.sum(window_1000[:, 1].astype(float) < 0) + np.sum(window_1000[:, 2].astype(float) < 0)

    window_1000[window_1000[:, 1].astype(float) < 0, 1] = 0
    window_1000[window_1000[:, 2].astype(float) < 0, 1] = 0

    size_of_matched_chr_vec = np.array([hg38[chrm] for chrm in window_1000[:, 0]], dtype = int)
    #if coordinate > chromosome size then change to chromosome size
    number_of_overbound = np.sum((window_1000[:, 1].astype(float) - size_of_matched_chr_vec) > 0) + np.sum((window_1000[:, 2].astype(float) - size_of_matched_chr_vec) > 0)

    window_1000[(window_1000[:, 1].astype(float) - size_of_matched_chr_vec) > 0, 1] = size_of_matched_chr_vec[(window_1000[:, 1].astype(float) - size_of_matched_chr_vec) > 0]
    window_1000[(window_1000[:, 2].astype(float) - size_of_matched_chr_vec) > 0, 2] = size_of_matched_chr_vec[(window_1000[:, 2].astype(float) - size_of_matched_chr_vec) > 0]

    print("Number of records with exceeding boundary: ", number_of_overbound)
    print("Number of records coordinate < 0:", lower_than_0)

    fname = input_bed.split("/")[-1]
    fname = fname.split(".")[0]
    cwd = os.getcwd()
    window_out_path = f"{out_dir}/{fname}_1000_window.bed"
    smaller_out_path = f"{out_dir}/{fname}_smaller_1000.bed"
    larger_out_path = f"{out_dir}/{fname}_larger_1000.bed"

    window_seq_name_list = seq_name_list = [":".join([bed[0], "-".join([bed[1], bed[2]])])for bed in window_1000]
    smaller_seq_name_list = seq_name_list = [":".join([bed[0], "-".join([bed[1], bed[2]])])for bed in smaller_1000]
    larger_seq_name_list = seq_name_list = [":".join([bed[0], "-".join([bed[1], bed[2]])])for bed in larger_1000]

    seq_name_bed(window_out_path, window_seq_name_list)
    seq_name_bed(smaller_out_path, smaller_seq_name_list)
    seq_name_bed(larger_out_path, larger_seq_name_list)

    print(f"Window written to: {window_out_path}")
    print(f"Smaller than 1000 written to: {smaller_out_path}")
    print(f"Larger than 1000 written to: {larger_out_path}")
    
def read_fasta(Dir):
    names = []
    seqs = []
    with open(Dir) as fp:
        for name, seq in get_fasta(fp):
            names.append(name[1:]) # remove >
            seqs.append(seq)
    return([names, seqs])

# convert sequence into np.array

def seq_to_array(seq_names, list_of_seqs, eccDNA = False):
    seq_list = []
    index_to_delete = []
    for index, seq in enumerate(list_of_seqs):
        one_hot = None
        one_hot = seq_to_one_hot(seq)
        if one_hot is None:
            index_to_delete.append(index)
            seq_list.append(np.zeros((1000, 4)))
            continue
        one_hot = np.array(one_hot, dtype = "int8")
        one_hot = np.reshape(one_hot, (1000, 4))
        seq_list.append(one_hot)
    seq_array = np.stack(seq_list)
    return(seq_names, seq_array, index_to_delete)

def seq_to_one_hot(seq, window_size = 1000):
    one_hot = pd.get_dummies(list(seq.upper()))
    if set("ACGT") != set(one_hot.columns): return None
    if one_hot.shape[0] == window_size: return(one_hot)
    #Append 0s to window_size
    appended = append_to_window(one_hot, window_size = window_size)
    return(appended)

def append_to_window(seq_array, window_size = 1000):
    #Append 0s to window_size
    appended = pd.DataFrame(seq_array, columns = ['A', 'C', 'G', 'T'])
    half = (window_size - (seq_array.shape[0]))/2 
    pre_num = math.floor(half)
    post_num = math.ceil(half)
    pre = pd.DataFrame([[0, 0, 0, 0]]*pre_num, columns = ['A', 'C', 'G', 'T'])
    post = pd.DataFrame([[0, 0, 0, 0]]*post_num, columns = ['A', 'C', 'G', 'T'])
    appended = pre.append(appended, ignore_index = True)
    appended = appended.append(post, ignore_index = True)
    return(appended)

def split_flanking_from_window(window_1000_array, eccDNA_array):
    flanking_array = window_1000_array - eccDNA_array
    result_arrray = np.dstack([eccDNA_array, flanking_array]) #Combine eccDNA and flanking (1000*4 to 1000*4)
    return(result_arrray)

def obtain_eccDNA_control(eccDNA_seqs_array ,control_1000_window_array):
    eccDNA_index_list = [np.apply_along_axis(sum, 1, seq) for seq in eccDNA_seqs_array]
    random.shuffle(eccDNA_index_list)
    control_eccDNA_list = []
    for index, window_1000 in enumerate(control_1000_window_array):
        try:
            control_eccDNA_seq_array = np.compress(eccDNA_index_list[index] ,window_1000, axis = 0)
        except IndexError:
            rand_num = random.sample(range(len(eccDNA_index_list)), 1)
            control_eccDNA_seq_array = np.compress(eccDNA_index_list[rand_num[0]] ,window_1000, axis = 0)
        control_eccDNA_list.append(np.array(append_to_window(control_eccDNA_seq_array), dtype = "int8")) #Append to 1000 window
    control_eccDNA_array = np.stack(control_eccDNA_list)
    return(control_eccDNA_array)

def eccDNA_window_identity(eccDNA_seqs_array, window_seqs_array, eccDNA_ind, window_ind):
    eccDNA = eccDNA_seqs_array[eccDNA_ind]
    eccDNA_length = sum(np.apply_along_axis(sum, 1, eccDNA))
    intersection = window_seqs_array[window_ind] == eccDNA
    intersection_length = sum(np.apply_along_axis(sum, 1, intersection) == 4)
    return(intersection_length == eccDNA_length)

def prob_to_discrete_prediction(predictions):
    predictions_01 = (predictions > 0.5).astype("int")
    list_prediction = []
    for i in predictions_01:
        if i[0] == 1:
            list_prediction.append(0)
            continue
        if i[1] == 1:
            list_prediction.append(1)
            continue
    return(list_prediction)

def test_and_confusion_mat(model, model_name, test_x, test_y_org, epi_type, cell, index = False):
    
    #model.load_weights('./checkpoints/'+model_name)
    model.compile(optimizer="adam",
        loss=tf.keras.losses.BinaryCrossentropy(),
        metrics=[
            tf.keras.losses.BinaryCrossentropy(from_logits=True, name='binary_crossentropy'),
            'accuracy'])
    predictions = model.predict(test_x)
    predictions_01 = (predictions > 0.5).astype("int")
    list_prediction = []
    for i in predictions_01:
        if i[0] == 1:
            list_prediction.append(0)
            continue
        if i[1] == 1:
            list_prediction.append(1)
            continue

    test_y_org = test_y_org.astype("int")
    pre_dict = {"Actual":test_y_org, "Prediction":list_prediction}
    pre_dict = pd.DataFrame(pre_dict)
    key = performance(list_prediction, test_y_org, epi_type, cell = cell)
    value = pd.crosstab(pre_dict.Prediction, pre_dict.Actual, margins = True)
    ind = [pre_dict.Prediction, pre_dict.Actual]
    if index == True:
        return(ind)
    else:
        return([key, value])

def confusion_mat(predictions, test_y_org, epi_type):
    test_y_org = test_y_org.astype("int")
    pre_dict = {"Actual":test_y_org, "Prediction":predictions}
    pre_dict = pd.DataFrame(pre_dict)
    key = performance(predictions, test_y_org, epi_type)
    value = pd.crosstab(pre_dict.Prediction, pre_dict.Actual, margins = True)
    return([key, value])

def performance(list_prediction, test_y_org, epi_type, cell):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for index, value in enumerate(list_prediction):
        if list_prediction[index] == test_y_org[index]:
            if list_prediction[index] == 1:
                TP += 1
            if list_prediction[index] == 0:
                TN += 1
        if list_prediction[index] != test_y_org[index]:
            if list_prediction[index] == 1:
                FP += 1
            if list_prediction[index] == 0:
                FN += 1
    sensitivity = round(TP/(TP+FN), 4)
    precision = round(TP/(TP+FP), 4)
    F1 = round(2/(1/sensitivity + 1/precision), 4)
    acc = round((TP+TN)/(TP+TN+FP+FN), 4)
    return(f"{cell}, {epi_type}, Sensitivity: {sensitivity}, Precision:{precision}, F1: {F1}, acc: {acc}")

def get_train_test(control_array, eccDNA_array, test_ratio = 0.2):
    train_num = round(min(len(control_array), len(eccDNA_array))*(1-test_ratio))
    test_num = round(min(len(control_array), len(eccDNA_array))*test_ratio)

    eccDNA_shuffle_ind = random.sample(range(len(eccDNA_array)), train_num)
    control_shuffle_ind = random.sample(range(len(control_array)), train_num)

    train_shuffle_ind = [eccDNA_shuffle_ind, control_shuffle_ind]
    
    train_x = np.concatenate([eccDNA_array[eccDNA_shuffle_ind], control_array[control_shuffle_ind]])
    train_x = train_x.astype(np.int8)
    
    train_eccDNA_y = [1]*train_num
    train_control_y = [0]*train_num
    train_y = train_eccDNA_y + train_control_y
    train_y = pd.get_dummies(train_y)
    train_y = np.array(train_y, dtype = "int8")
    
    test_eccDNA_shuffle_ind = random.sample(set(range(len(eccDNA_array))) - set(eccDNA_shuffle_ind), test_num)
    test_control_shuffle_ind = random.sample(set(range(len(control_array))) - set(control_shuffle_ind), test_num) 
    
    test_shuffle_ind = [test_eccDNA_shuffle_ind, test_control_shuffle_ind]
    
    test_x = np.concatenate([eccDNA_array[test_eccDNA_shuffle_ind], control_array[test_control_shuffle_ind]])
    test_x = test_x.astype(np.int8)

    test_eccDNA_y = [1]*test_num
    test_control_y = [0]*test_num
    test_y = test_eccDNA_y + test_control_y
    test_y = pd.get_dummies(test_y)
    test_y = np.array(test_y, dtype = "int8")
    test_y_org = np.concatenate([[1]*test_num, [0]*test_num],)
    test_y_org = np.array(test_y_org, dtype = "int8")
    return(train_x, train_y, test_x, test_y, test_y_org, train_shuffle_ind, test_shuffle_ind)

def get_callbacks(name, out_dir):
    cwd = os.getcwd()
    return [
        tf.keras.callbacks.ModelCheckpoint(
        filepath = f"{out_dir}/models/"+name,
        monitor="val_binary_crossentropy",
        verbose=1,
        save_best_only=True,
        save_weights_only=True,
        mode="min",
        save_freq="epoch"
        )
  ]

def get_optimizer():
    return tf.keras.optimizers.Adam(learning_rate=0.001)

def compile_and_fit(name, train_x, train_y, out_dir, optimizer=None, epochs=100, batch_size = 32):
    if optimizer is None:
        optimizer = get_optimizer()
    n_timesteps, n_features, n_outputs = train_x.shape[1], train_x.shape[2], 2
    model = Sequential([
        Conv1D(filters=16, kernel_size=3, activation='relu', kernel_regularizer=regularizers.l2(0.001) ,input_shape=(n_timesteps,n_features)),
        Dropout(0.2),
        Conv1D(filters=16, kernel_size=3, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
        Dropout(0.2),
        MaxPool1D(pool_size=2),
        Flatten(),
        Dense(100, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
        Dropout(0.2),
        Dense(n_outputs, activation = 'softmax')
    ])
    model.compile(optimizer=optimizer,
        loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
        metrics=[
            tf.keras.losses.BinaryCrossentropy(
            from_logits=True, name='binary_crossentropy'),
            'accuracy'])
    model.summary()
    model.fit(
        train_x,
        train_y,
        epochs = epochs,
        validation_split = 0.2,
        batch_size = batch_size,
        callbacks = get_callbacks(name, out_dir = out_dir),
        verbose = 0)
    return(model)

def split_data_to_chunks(seq_name_pairs, num_chunks = 8):
    if len(seq_name_pairs[0]) != len(seq_name_pairs[0]):
        print("pairs with different length!")
        return
    length = len(seq_name_pairs[0])
    size_chunks = len(seq_name_pairs[0])//num_chunks
    names = seq_name_pairs[0]
    seqs = seq_name_pairs[1]
    if length % size_chunks == 0:
        splitted = [tuple([names[start*size_chunks: min((start + 1)*size_chunks, length)],
                          seqs[start*size_chunks: min((start + 1)*size_chunks, length)]])
                          for start in range(num_chunks)]
    else:
        splitted = [tuple([names[start*size_chunks: min((start + 1)*size_chunks, length)],
                          seqs[start*size_chunks: min((start + 1)*size_chunks, length)]])
                          for start in range(num_chunks+1)]
    return(splitted)

def delete_name_seq_pair_with_ind(names, seq_array, del_index):
    remained_names = np.delete(np.array(names), del_index, axis = 0)
    remained_seqs = np.delete(seq_array, del_index, axis = 0)
    return(remained_names, remained_seqs)
def delete_seq_with_N(eccDNA_seqs_names, eccDNA_seqs_array, window_seqs_names, window_seqs_array, window_del_ind):
    return(delete_name_seq_pair_with_ind(eccDNA_seqs_names, eccDNA_seqs_array, window_del_ind),
           delete_name_seq_pair_with_ind(window_seqs_names, window_seqs_array, window_del_ind))

#read eccDNA from fa file and convert them into tensors
def split_array(eccDNA_name_seq_array_slice, window_name_seq_array_slice, control_name_seq_array_slice):
    eccDNA_seqs_names, eccDNA_seqs_array, eccDNA_del_ind = eccDNA_name_seq_array_slice
    window_seqs_names, window_seqs_array, window_del_ind = window_name_seq_array_slice
    control_window_names, control_window_array, control_del_ind = control_name_seq_array_slice

    eccDNA_name_seq_pair, window_name_seq_pair = delete_seq_with_N(eccDNA_seqs_names, eccDNA_seqs_array, window_seqs_names, window_seqs_array, window_del_ind)
    
    control_window_names = np.delete(np.array(control_name_seq_array_slice[0]), control_del_ind, axis = 0)
    control_window_array = np.delete(control_name_seq_array_slice[1], control_del_ind, axis = 0)
    control_eccDNA_array = obtain_eccDNA_control(eccDNA_name_seq_pair[1] ,control_window_array)
    
    control_flank_split_array = split_flanking_from_window(control_window_array ,control_eccDNA_array)
    eccDNA_flank_split_array = split_flanking_from_window(window_name_seq_pair[1], eccDNA_name_seq_pair[1])
    return(control_flank_split_array, eccDNA_flank_split_array, control_window_names, window_name_seq_pair[0])

def save_prediction_result(fname, test_seq_names, prediction_result):
    with open(fname, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["Genomic coordinates", "Probability (non-eccDNA)", "Probability (eccDNA)"])
        for ind, prob in enumerate(prediction_result):
            writer.writerow([test_seq_names[ind], prob[0], prob[1]])
    print(f"Results written to: {fname}")

def gen_pred_pos_neg(fname, test_seq_names, prediction_result, genome_fa):
    pos_bed = f"{fname}_pred_positive.bed"
    neg_bed = f"{fname}_pred_negative.bed"
    pos_fa = f"{fname}_pred_positive.fa"
    neg_fa = f"{fname}_pred_negative.fa"
    with open(pos_bed, 'w', newline='') as pos_out:
        with open(neg_bed, 'w', newline='') as neg_out:
            writer_pos_out = csv.writer(pos_out, delimiter='\t')
            writer_neg_out = csv.writer(neg_out, delimiter='\t')
            for ind, prob in enumerate(prediction_result):
                chrm = test_seq_names[ind].split(":")[0]
                start, end = test_seq_names[ind].split(":")[1].split("-")
                seq_name_sep = [chrm, start, end]
                if prob[1] >= 0.5:
                    writer_pos_out.writerow(seq_name_sep)
                else:
                    writer_neg_out.writerow(seq_name_sep)
    #Convert bed files to fasta files
    os.system(f"bedtools getfasta -fi {genome_fa} -bed {pos_bed} -fo {pos_fa}")
    os.system(f"bedtools getfasta -fi {genome_fa} -bed {neg_bed} -fo {neg_fa}")

    print(f"Results written to: {pos_bed}, {neg_bed}, {pos_fa}, {neg_fa}")
            
def seq_name_bed(fname, seq_name_list):
    
    with open(fname, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        for seq_name in seq_name_list:
            chrm = seq_name.split(":")[0]
            start, end = seq_name.split(":")[1].split("-")
            seq_name_sep = [chrm, start, end]
            writer.writerow(seq_name_sep)

#Preprocess genome fasta, eccDNA bed files to generate positive and negative label data
def preprocessing(ncpus, genome_fa, genome, genome_gap, eccdna_bed_dir, out_dir):
    
    #Convert eccDNA bed to 1000 bp window bed
    fname = eccdna_bed_dir.split("/")[-1]
    fname = fname.split(".")[0]
    cwd = os.getcwd()

    #Check if output directory exists, if not, make dir
    if not os.path.isdir(out_dir):
        os.system(f"mkdir {out_dir}")
    get_1000_window(eccdna_bed_dir, genome, out_dir)
    #Check if there's output directory
    #Shuffle the 1000 bp window
    window_bed_Dir = f'{out_dir}/{fname}_1000_window.bed'
    Con_bed_Dir = f'{out_dir}/rand_{fname}_1000_window.bed'
    small_eccDNA_bed_Dir = f'{out_dir}/{fname}_smaller_1000.bed'
    excl_bed_Dir = f'{out_dir}/{fname}_1000_window_excl.bed'
    os.system(f"cat {genome_gap} {window_bed_Dir} | sort -k1,1 -k2,2n > {excl_bed_Dir}")
    os.system(f"bedtools shuffle -i {window_bed_Dir} -g {genome} -excl {excl_bed_Dir} -maxTries 5000  > {Con_bed_Dir}")

    #Convert bed to fasta
    bed_list = [window_bed_Dir, small_eccDNA_bed_Dir, Con_bed_Dir]
    fa_list = [bed.split(".bed")[0] + ".fa" for bed in bed_list]
    for i, fa in enumerate(fa_list):
        os.system(f"bedtools getfasta -fi {genome_fa} -bed {bed_list[i]} -fo {fa}")

    #Convert fasta to numpy array
    window_names, window_seqs = read_fasta(fa_list[0])
    eccDNA_names_smaller_1000, eccDNA_seqs_smaller_1000 = read_fasta(fa_list[1])
    control_names, control_1000_window_seqs_list = read_fasta(fa_list[2])
    start = time.time()
    eccDNA_inputs = split_data_to_chunks([eccDNA_names_smaller_1000, eccDNA_seqs_smaller_1000], num_chunks = ncpus)
    window_inputs = split_data_to_chunks([window_names, window_seqs], num_chunks = ncpus)
    control_inputs = split_data_to_chunks([control_names, control_1000_window_seqs_list], num_chunks = ncpus)
    pool = Pool(ncpus)
    eccDNA_name_seq_array = pool.starmap(seq_to_array, eccDNA_inputs)
    window_name_seq_array = pool.starmap(seq_to_array, window_inputs)
    control_name_seq_array = pool.starmap(seq_to_array, control_inputs)
    eccDNA_window_control_combined = [tuple([eccDNA_name_seq_array[i], window_name_seq_array[i], 
                                        control_name_seq_array[i]]) for i in range(len(eccDNA_name_seq_array))]
    flank_split_array = pool.starmap(split_array, eccDNA_window_control_combined)

    control_flank_split_array = np.concatenate([flank_split_array[i][0] for i in range(len(flank_split_array))])
    eccDNA_flank_split_array = np.concatenate([flank_split_array[i][1] for i in range(len(flank_split_array))])
    control_window_names = np.concatenate([flank_split_array[i][2] for i in range(len(flank_split_array))])
    window_seqs_names = np.concatenate([flank_split_array[i][3] for i in range(len(flank_split_array))])

    control_window_names, uniq_control_ind = np.unique(control_window_names, return_index=True)
    window_seqs_names, uniq_window_ind = np.unique(window_seqs_names, return_index=True)
    control_flank_split_array = control_flank_split_array[uniq_control_ind]
    eccDNA_flank_split_array = eccDNA_flank_split_array[uniq_window_ind]
    end = time.time()
    elapsed_time = (end -start)/60
    print(f"Finished converting {fname} DNA seq to tensors, took {elapsed_time} min")
    return([control_flank_split_array, eccDNA_flank_split_array, control_window_names, window_seqs_names])

def fa_to_arr(ncpus, window_fa_Dir, small_eccDNA_fa_Dir, Con_fa_Dir):

    #Convert fasta to numpy array
    window_names, window_seqs = read_fasta(window_fa_Dir)
    eccDNA_names_smaller_1000, eccDNA_seqs_smaller_1000 = read_fasta(small_eccDNA_fa_Dir)
    control_names, control_1000_window_seqs_list = read_fasta(Con_fa_Dir)
    start = time.time()
    eccDNA_inputs = split_data_to_chunks([eccDNA_names_smaller_1000, eccDNA_seqs_smaller_1000], num_chunks = ncpus)
    window_inputs = split_data_to_chunks([window_names, window_seqs], num_chunks = ncpus)
    control_inputs = split_data_to_chunks([control_names, control_1000_window_seqs_list], num_chunks = ncpus)
    pool = Pool(ncpus)
    eccDNA_name_seq_array = pool.starmap(seq_to_array, eccDNA_inputs)
    window_name_seq_array = pool.starmap(seq_to_array, window_inputs)
    control_name_seq_array = pool.starmap(seq_to_array, control_inputs)
    eccDNA_window_control_combined = [tuple([eccDNA_name_seq_array[i], window_name_seq_array[i], 
                                        control_name_seq_array[i]]) for i in range(len(eccDNA_name_seq_array))]
    flank_split_array = pool.starmap(split_array, eccDNA_window_control_combined)

    control_flank_split_array = np.concatenate([flank_split_array[i][0] for i in range(len(flank_split_array))])
    eccDNA_flank_split_array = np.concatenate([flank_split_array[i][1] for i in range(len(flank_split_array))])
    control_window_names = np.concatenate([flank_split_array[i][2] for i in range(len(flank_split_array))])
    window_seqs_names = np.concatenate([flank_split_array[i][3] for i in range(len(flank_split_array))])

    control_window_names, uniq_control_ind = np.unique(control_window_names, return_index=True)
    window_seqs_names, uniq_window_ind = np.unique(window_seqs_names, return_index=True)
    control_flank_split_array = control_flank_split_array[uniq_control_ind]
    eccDNA_flank_split_array = eccDNA_flank_split_array[uniq_window_ind]
    end = time.time()
    elapsed_time = (end -start)/60
    print(f"Finished converting {fname} DNA seq to tensors, took {elapsed_time} min")
    return([control_flank_split_array, eccDNA_flank_split_array, control_window_names, window_seqs_names])

def eccDNA_CNN_train(control_flank_split_array,
                     eccDNA_flank_split_array,
                     control_window_names,
                     window_seqs_names,
                     eccdna_bed_dir,
                     test_ratio,
                     epochs,
                     batch_size,
                     output_bed):
    ##For model training
    #Split training and testing data
    train_x, train_y, test_x, test_y, test_y_org, train_shuffle_ind, test_shuffle_ind = get_train_test(control_flank_split_array, eccDNA_flank_split_array, test_ratio = test_ratio)

    #Train a model
    train_x_window = eccDNA_flank_split_array[train_shuffle_ind[0]]
    train_x_control = control_flank_split_array[train_shuffle_ind[1]]    
    train_x, train_y, _, _, _, _, _ = get_train_test(train_x_control, train_x_window, test_ratio = 0)

    fname = eccdna_bed_dir.split("/")[-1]
    fname = fname.split(".")[0]
    model_name = fname + "_CNN_model"
    model = compile_and_fit(model_name, train_x, train_y, out_dir = output_bed,optimizer=None, epochs = int(epochs), batch_size = int(batch_size))
    
    #Output confusion matrix
    confusion_matrix = test_and_confusion_mat(model, model_name, test_x, test_y_org, epi_type= "testing result", cell = model_name)
    print(confusion_matrix[0])
    print(confusion_matrix[1])
    cwd = os.getcwd()
    
    if isinstance(output_bed, str):
        #Output bed for training and testing windows
        train_positive = window_seqs_names[train_shuffle_ind[0]]
        train_negative = control_window_names[train_shuffle_ind[1]]
        test_positive = window_seqs_names[test_shuffle_ind[0]]
        test_negative = control_window_names[test_shuffle_ind[1]]
        seq_name_bed(f"{output_bed}/{fname}_eccdna_train_positive.bed", train_positive)
        seq_name_bed(f"{output_bed}/{fname}_eccdna_train_negative.bed", train_negative)
        seq_name_bed(f"{output_bed}/{fname}_eccdna_test_positive.bed", test_positive)
        seq_name_bed(f"{output_bed}/{fname}_eccdna_test_negative.bed", test_negative)

##For prediction based on pre-trained models
def eccDNA_CNN_predict(control_flank_split_array,
                       eccDNA_flank_split_array,
                       control_window_names,
                       window_seqs_names,
                       eccdna_bed_dir,
                       chosen_model,
                       genome_fa,
                       output_prob):
    fname = eccdna_bed_dir.split("/")[-1]
    fname = fname.split(".")[0]
    
    _, _, test_x, test_y, test_y_org, _, test_shuffle_ind = get_train_test(control_flank_split_array, eccDNA_flank_split_array, test_ratio = 1)
    n_channels = 8
    model = Sequential([
            Conv1D(filters=16, kernel_size=3, activation='relu', kernel_regularizer=regularizers.l2(0.001) ,input_shape=(1000, n_channels)),
            Dropout(0.2),
            Conv1D(filters=16, kernel_size=3, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            Dropout(0.2),
            MaxPool1D(pool_size=2),
            Flatten(),
            Dense(100, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            Dropout(0.2),
            Dense(2, activation = 'softmax')
        ])
    cwd = os.getcwd()
    #Load pre-trained model
    model_name = f"{chosen_model}_seq_only_ensemble_base_tf115"
    model.load_weights(f"{cwd}/pre_trained_models/{model_name}")
    prediction = model.predict(test_x)
    prediction_window = prediction[:int(prediction.shape[0]/2)]

    if isinstance(output_prob, str):
        #Write prediction bed files
        test_seq_names = window_seqs_names[test_shuffle_ind[0]]
        save_prediction_result(f"{output_prob}/{chosen_model}/{fname}_prob.tsv", test_seq_names, prediction_window)
        gen_pred_pos_neg(f"{output_prob}/{chosen_model}/{fname}", test_seq_names, prediction_window, genome_fa)
        
def gen_pos_neg_seq(pred_label_dir, pred_seq_dir):
    pos_ind_list = []
    neg_ind_list = []

    #Parsing prediction label
    try:
        with open(pred_label_dir) as in_f:
            reader = csv.reader(in_f, delimiter = "\t")
            for ind, label in reader:
                if len(ind) == 0:
                    continue
                if label == "TP" or label == "FP":
                    pos_ind_list.append(ind)
                if label == "TN" or label == "FN":
                    neg_ind_list.append(ind)
    except FileNotFoundError:
        print("Please predict eccDNAs first")

    #Parsing prediction sequences and writing positive and negative results
    out_pos_dir = pred_seq_dir.split(".tsv")[0] + "_pred_positive.fa"
    out_neg_dir = pred_seq_dir.split(".tsv")[0] + "_pred_negative.fa"
    with open(out_pos_dir, "w") as out_pos_seq:
        with open(out_neg_dir, "w") as out_neg_seq:
            pos_writer = csv.writer(out_pos_seq)
            neg_writer = csv.writer(out_neg_seq)
            with open(pred_seq_dir) as in_seq:
                reader = csv.reader(in_seq, delimiter = "\t")
                for ind, seq in reader:
                    if len(ind) == 0:
                        continue
                    if ind in pos_ind_list:
                        pos_writer.writerow([f">{ind}"])
                        pos_writer.writerow([seq])
                    if ind in neg_ind_list:
                        neg_writer.writerow([f">{ind}"])
                        neg_writer.writerow([seq])
    print(f"Written to: {out_pos_dir}, {out_neg_dir}")
    return([out_pos_dir, out_neg_dir])
        