import matplotlib.pyplot as plt
import os 
import sys
import numpy as np
import pandas as pd

from calculate_modularity import calculate_modularity
from null_model_und_sign import null_model_und_sign
from construct_nulls import construct_NG_null
from get_similarity_consensus import get_similarity_consensus

from FISHnet_support_functions import convert_tuples
from FISHnet_support_functions import smoove
from FISHnet_support_functions import Modularity_maximization
from FISHnet_support_functions import find_small_groups
from FISHnet_support_functions import adjust_small_groups
from FISHnet_support_functions import Size_exclude_communities
from FISHnet_support_functions import count_switches
from FISHnet_support_functions import find_close_values
from FISHnet_support_functions import find_index
from FISHnet_support_functions import find_plateau_points
from FISHnet_support_functions import acc_check 


def FISHnet_main(input_matrix,distance = [250],plateau_size = 4,window_size = 2,size_exclusion = 3, merge = 3):
    
    if plateau_size == 0 and len(distance) != 1:
        raise Exception("When having plateau_size = 0, please have only one distance in the list")
    if plateau_size != 0 and len(distance) == 1:
        raise Exception("When having plateau_size = 0, please have only one distance in the list")
   

    
    # Get information abount the input matrix
    walks = np.shape(input_matrix)[0]
    
    

    #Step 1 Conducting the threshold platue sweep
    ########################################################################################################
    avg_communities = []
    
    Max_size_of_matrix = np.shape(input_matrix)[1]

    max_search = int(np.nanmax(input_matrix))

    thresh_to_use = []

    for thresh in distance:
        #print(thresh)

        matrix2 = input_matrix < thresh

        matrix = smoove(matrix2,window=window_size)

        P = construct_NG_null(matrix)


        num_com = []
        for itr in range(20):
            q,communities = calculate_modularity(matrix,P,1,seed=None)
            communities_itr = len(np.unique(communities))


            ##########################
            if communities_itr == 2:
                communities_switches = count_switches(communities)
                if communities_switches > 2:
                    continue
                else:
                    communities_itr = communities_switches
            ##########################
            num_com.append(communities_itr)
        
        if len(num_com) == 0:
            mean_v = np.nan 
        else:
            mean_v = np.mean(num_com)
                
        rounded_v = np.round(mean_v,0)
        if rounded_v == Max_size_of_matrix:
            continue
        avg_communities.append(rounded_v)
        thresh_to_use.append(thresh)
    ########################################################################################################


    
    #Step 2 Get the thresholds we are going to use
    ########################################################################################################
    
    gammas_to_use, comm_to_have,group_dic = find_plateau_points(thresh_to_use,avg_communities,platue_len=plateau_size)
    
    ########################################################################################################
    
    
    
    #Step 3 Consense Calls from the same Platue Community
    ########################################################################################################
    Number_groups = len(group_dic.keys())

    Final_communities = np.zeros((Number_groups,Max_size_of_matrix))
    
    Per_plat_community = {}

    ng_counter = 0
    for item in group_dic.keys():
        consensuss = np.zeros((len(group_dic[item]),Max_size_of_matrix))
        counter = 0

        for thresh in group_dic[item]:
            #print(thresh)
            other = Modularity_maximization(input_matrix,threshold=thresh,window_size = window_size)
            consensuss[counter,:] = other
            counter+=1
        
        Per_plat_community[item] = consensuss # NEW

        fc,other = get_similarity_consensus(consensuss)

        Final_communities[ng_counter,:] = fc
        ng_counter+=1
       
    ########################################################################################################
    
    
    
    #Step 4 Size Exclusion (remove communities that are smaller than a specified value)
    ########################################################################################################
    
    Final_communities_new = np.zeros(np.shape(Final_communities))

    for itr in range(np.shape(Final_communities)[0]):
        fc = Final_communities[itr,:]
        New_corrected = Size_exclude_communities(fc,threshold_value=size_exclusion)
        Final_communities_new[itr,:] = New_corrected
        

    ########################################################################################################



    #Step 5 Replace Final communities with the size removed communities
    ########################################################################################################
    
    Final_communities = Final_communities_new.copy()
   
    ########################################################################################################
    
    
    
    #Step 6 Identify Boundaries from the Final communities 
    ########################################################################################################
    
    tot_new = {}
    for itr in range(np.shape(Final_communities)[0]):

        boundaries = {0:100}
        initial = Final_communities[itr,:][0]
        #print( Final_communities[itr,:])
        counter = 0
        for item in Final_communities[itr,:]:
            if item != initial:
                size = sum(Final_communities[itr,:] == initial)
                if size >= 0:
                    boundaries[counter] = size
                initial = item
            counter+=1

        boundaries[walks]=100

        tot_new[itr] = list(boundaries.keys())

    
    ########################################################################################################


    #Step 7 MERGE Boundaries Close to eachother!
    ########################################################################################################
    
    combined_boundaries = []
    for keys,items in tot_new.items():
        for values in items:
            if (values == 0) or (values==Max_size_of_matrix):
                continue
            if values not in combined_boundaries:
                combined_boundaries.append(values)


    result = find_close_values(combined_boundaries, threshold=merge)

    

    final_dic = {}

    for keys,items in tot_new.items():
        nw = []
        for item in items:

            if (item == 0) or (item == Max_size_of_matrix):
                nw.append(item)
                continue

            index = find_index(item,result)
            mean_value = np.mean(result[index])

            nw.append(mean_value)

        final_dic[keys] = nw
                
    
    ########################################################################################################
    
    
    # Convert Boundaries back to Domains
    
    Domain_output = convert_tuples(final_dic)
    
    # Output Dictionary of Final Domain Calls
    
    return Domain_output,group_dic

