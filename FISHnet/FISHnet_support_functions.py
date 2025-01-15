import matplotlib.pyplot as plt
import os 
import sys
import numpy as np
import pandas as pd

from calculate_modularity import calculate_modularity
from null_model_und_sign import null_model_und_sign
from construct_nulls import construct_NG_null
from get_similarity_consensus import get_similarity_consensus


def acc_check(call,Array,border = 2):
    
    new_array = []
    
    for item in Array:
        
        min_range = item - border
        
        max_range = item + border 
        
        if min_range <= call <= max_range:
            
            return 1
        
    return 0

def smoove(matrix,window):
    
    if window == 0:
        return matrix
    
    walks = np.shape(matrix)[0]
    
    Combined_bin = np.zeros((walks,walks))
    for i in range(walks):
        for j in range(walks):

            lf = i-window
            uf = i+window
            ls = j-window
            us = j+window

            if lf <0:
                lf = 0
            if ls <0:
                ls = 0
            if uf > walks:
                uf = walks
            if us > walks:
                us = walks
  
            mean_in = np.nanmean(matrix[lf:uf,ls:us])
            Combined_bin[i][j] = mean_in
            
    return Combined_bin



def Modularity_maximization(input_matrix,threshold=250,window_size = 2):
    walks = np.shape(input_matrix)[0]
    matrix2 = input_matrix < threshold
    
    matrix = smoove(matrix2,window=window_size)
    
    P = construct_NG_null(matrix) 
    
    ##########################################################
   
    gammas_to_use = [1]
    
    shapo = np.shape(matrix)[0]
    Final_communities = np.zeros((len(gammas_to_use),shapo))
    counter = 0
    for gammy in gammas_to_use:
        consensus = np.zeros((20,shapo))

        for itr in range(20):
            q,communities = calculate_modularity(matrix,P,gammy,seed=None)
            consensus[itr,:] = communities

        fc,other = get_similarity_consensus(consensus)
        Final_communities[counter,:] = fc
        counter+=1
        
    ##########################################################
    
            
    return Final_communities 



def find_small_groups(arr, threshold):
    small_groups = []
    current_group = []
    for i in range(len(arr)):
        if arr[i] == 1:
            current_group.append(i)
        else:
            if len(current_group) > 0:
                if len(current_group) <= threshold:
                    small_groups.append(current_group)
                current_group = []
    
    # Check if the last group is small
    if len(current_group) > 0 and len(current_group) <= threshold:
        small_groups.append(current_group)
    
    return small_groups





def adjust_small_groups(fc,groups):
    
    min_boundary = 0
    
    max_boundary = np.shape(fc)[0] - 1

    adjusted_groups = []

    for i in groups:
        mid_point = int(len(i)/2)
        counter = 0
        holder = []
        
        
        if min_boundary in i:
            adjusted_groups.append(list(np.array(i)+1))
            
            continue
        
        if max_boundary in i:
            adjusted_groups.append(list(np.array(i)-1))
            continue
        
        for item in i:
                  
            if counter < mid_point:
                holder.append(item - 1)

            else:
                holder.append(item+1)

            counter+=1

        adjusted_groups.append(holder)
        
    return adjusted_groups




def Size_exclude_communities(community,threshold_value):
    
    fc = community.copy()
    max_size = np.shape(fc)[0]
    for item in np.unique(fc):
        
        
        
        binary = (fc == item)*1

        groups = find_small_groups(binary,threshold=threshold_value)

        adjusted_groups = adjust_small_groups(fc,groups)

        for item,aItem in zip(groups,adjusted_groups):
            
            

            if max_size in item:
                for bound,abound in zip(item,aItem):
                    fc[bound] = fc[abound]
                continue


            if 0 in item:
                for bound,abound in zip(reversed(item),reversed(aItem)):
                     fc[bound] = fc[abound]
                continue


            ############################################################################
            N = int(np.round(len(item)/2))

            first = 0
            last = N
            for indexer in range(N):
                bound = item[first]
                bound_end = item[last]

                abound = aItem[first]
                abound_end = aItem[last]

                fc[bound] = fc[abound]
                fc[bound_end] = fc[abound_end]

                first+=1
                last-=1
            ############################################################################
            
    return fc



def count_switches(arr):
    switch_count = 1
    prev_num = None

    for num in arr:
        if prev_num is not None and num != prev_num:
            switch_count += 1
        prev_num = num

    return switch_count




def find_close_values(lst, threshold=2):
    
    if len(lst) == 0:
        return []
    
    lst.sort()  # Sort the list to simplify the comparison process
    groups = []
    
    current_group = [lst[0]]

    for i in range(1, len(lst)):
        if abs(lst[i] - lst[i-1]) <= threshold:
            current_group.append(lst[i])
        else:
            groups.append(current_group)
            current_group = [lst[i]]

    groups.append(current_group)  # Append the last group

    return groups


def find_index(value,list_of_list):
    indices = None
    counter = 0
    for item in list_of_list:
        if value in item:
            indices = counter
        counter+=1
        
    return indices



def find_plateau_points(datax,datay,platue_len = 4):
    
    # case when plat is 0
    if platue_len == 0:
        return datax,datay,{0:datax}
    
    # case when all values are same community num
    #v1 = datay[0]
    #are_all_same = np.sum(np.array(datay) == v1)
    
    #if (are_all_same == len(datay)) and (len(datay) >= platue_len):
    #    return datax,datay,{0:datax}
        
    
    
    plateau_points_x = []
    plateau_points_y = []
    group_dic = {}
    plateau_points_dic = []
    
    current_x = None
    current_y = datay[0]
    plateau_length = 0
    
    builder_plats_x = []
    builder_plats_y = []
    counter =0
    for x, y in zip(datax,datay):
        
        
        if y == current_y:
            plateau_length += 1
            builder_plats_x.append(x)
            builder_plats_y.append(y)
            
            
        else:
            if (plateau_length >= platue_len):
                group_dic[counter] = plateau_points_dic
                plateau_points_dic = []
                counter+=1
            plateau_length = 1
            current_x = x
            current_y = y
            builder_plats_x = [current_x]
            builder_plats_y = [current_y]
            

        if plateau_length >= platue_len:
            plateau_points_x = plateau_points_x + builder_plats_x
            plateau_points_y = plateau_points_y + builder_plats_y
            plateau_points_dic = plateau_points_dic + builder_plats_x
            builder_plats_x = []
            builder_plats_y = []
            group_dic[counter] = plateau_points_dic
    
    return plateau_points_x,plateau_points_y,group_dic


def convert_tuples(domain_dict):
    filtered_domains = {}

    for key, boundaries in domain_dict.items():
        # Initialize a list to store valid domain pairs for the current key
        valid_domains = []
        
        # Create domain pairs and filter them based on the min_size
        for i in range(len(boundaries) - 1):
            start, end = boundaries[i], boundaries[i + 1]
            if end - start >= 0:
                valid_domains.append((start, end))
        
        # Update the dictionary with the filtered domains
        filtered_domains[key] = valid_domains

    return filtered_domains



