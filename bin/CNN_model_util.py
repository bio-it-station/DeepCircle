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

def get_callbacks(name):
  return [
    tf.keras.callbacks.ModelCheckpoint(
    filepath = "../output/models/"+name,
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

def compile_and_fit(name, train_x, train_y, optimizer=None, epochs=100, batch_size = 32):
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
        callbacks = get_callbacks(name),
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
            
def seq_name_bed(fname, seq_name_list):
    with open(fname, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        for seq_name in seq_name_list:
            chrm = seq_name.split(":")[0]
            start, end = seq_name.split(":")[1].split("-")
            seq_name_sep = [chrm, start, end]
            writer.writerow(seq_name_sep)

#Preprocessing
def preprocessing(ncpus, genome_fa, genome, eccdna_bed_dir):
    #Convert eccDNA bed to 1000 bp window bed
    
    match = re.search("(\w+\.bed)", eccdna_bed_dir)
    fname = match.group()
    fname = fname.split(".")[0]
    os.system(f"Rscript get_1000_window.R {eccdna_bed_dir} {genome}")

    #Shuffle the 1000 bp window
    window_bed_Dir = f'../output/{fname}_1000_window.bed'
    Con_bed_Dir = f'../output/rand_{fname}_1000_window.bed'
    small_eccDNA_bed_Dir = f'../output/{fname}_smaller_1000.bed'
    os.system(f"bedtools shuffle -i {window_bed_Dir} -g {genome} -excl {window_bed_Dir} -maxTries 5000  > {Con_bed_Dir}")

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

def eccDNA_CNN_train(control_flank_split_array,
                     eccDNA_flank_split_array,
                     control_window_names,
                     window_seqs_names,
                     eccdna_bed_dir,
                     test_ratio,
                     epochs,
                     batch_size,
                     output_bed = True):
    ##For model training
    #Split training and testing data
    train_x, train_y, test_x, test_y, test_y_org, train_shuffle_ind, test_shuffle_ind = get_train_test(control_flank_split_array, eccDNA_flank_split_array, test_ratio = test_ratio)

    #Train a model
    train_x_window = eccDNA_flank_split_array[train_shuffle_ind[0]]
    train_x_control = control_flank_split_array[train_shuffle_ind[1]]    
    train_x, train_y, _, _, _, _, _ = get_train_test(train_x_control, train_x_window, test_ratio = 0)
    match = re.search("(\w+\.bed)", eccdna_bed_dir)
    fname = match.group()
    fname = fname.split(".")[0]
    model_name = fname + "_CNN_model"
    model = compile_and_fit(model_name, train_x, train_y, optimizer=None, epochs = int(epochs), batch_size = int(batch_size))
    
    #Output confusion matrix
    confusion_matrix = test_and_confusion_mat(model, model_name, test_x, test_y_org, epi_type= "testing result", cell = model_name)
    print(confusion_matrix[0])
    print(confusion_matrix[1])

    if output_bed:
        #Output bed for training and testing windows
        train_positive = window_seqs_names[train_shuffle_ind[0]]
        train_negative = control_window_names[train_shuffle_ind[1]]
        test_positive = window_seqs_names[test_shuffle_ind[0]]
        test_negative = control_window_names[test_shuffle_ind[1]]
        seq_name_bed(f"../output/{fname}_eccdna_train_positive.bed", train_positive)
        seq_name_bed(f"../output/{fname}_eccdna_train_negative.bed", train_negative)
        seq_name_bed(f"../output/{fname}_eccdna_test_positive.bed", test_positive)
        seq_name_bed(f"../output/{fname}_eccdna_test_negative.bed", test_negative)

##For prediction based on pre-trained models
def eccDNA_CNN_predict(control_flank_split_array,
                       eccDNA_flank_split_array,
                       control_window_names,
                       window_seqs_names,
                       eccdna_bed_dir,
                       chosen_model,
                       output_prob):
    match = re.search("(\w+\.bed)", eccdna_bed_dir)
    fname = match.group()
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
    
    #Load pre-trained model
    model_name = f"{chosen_model}_seq_only_ensemble_base_tf115"
    model.load_weights(f"../pre_trained_models/{model_name}")
    prediction = model.predict(test_x)
    prediction_window = prediction[:int(prediction.shape[0]/2)]
    if output_prob:
        #Write prediction bed files
        test_seq_names = window_seqs_names[test_shuffle_ind[0]]
        save_prediction_result(f"../output/{fname}_prob.tsv", test_seq_names, prediction_window)
        