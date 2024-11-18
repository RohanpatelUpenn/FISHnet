Howdy! Within this folder you will find the FISHnet tutorial. FISHnet is a graph based algorithm to find chromatin domains within sequential oligopaints data. See https://www.biorxiv.org/content/10.1101/2024.06.18.599627v1.article-info for more details on how FISHnet works. 

The run time for this tutorial is ~5 minutes. 

Within this folder you will find:


1) Python scripts that are used to run FISHnet (housed in the folder called `FISHnet`). Scripts include the following:   
    - FISHnet_main.py.
    - FISHnet_support_functions.py
    - genlouvain.py
    - construct_nulls.py
    - randmio_und_signed.py
    - calculate_modularity.py
    - get_similarity_consensus.py
    - null_model_und_sign.py

2) The FISHnet_Tutorial_notebook.ipynb in a juypter notebook.

3) HCT116_chr21_X.txt which is the Bintu 2018 HCT116 Chr21:34.6-37.1 Mb Dataset used in the tutorial.


Make sure you have the proper python packages installed as dictated by requirments.txt before running the tutorial! Have fun and happy domain calling!
