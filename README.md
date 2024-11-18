# FISHnet - a graph theory approach for calling chromatin domains in sequential Oligopaints data

![FISHnet_video_twitter](https://github.com/user-attachments/assets/e5f2aff0-ad25-4568-83d6-3082c8dc7f60)

## What is FISHnet?:

FISHnet is graph theory method to detect chromatin domains and boundaries in pairwise distance matrices from sequential Oligopaints data. FISHnet uses various distance thresholds as a resolution parameter to find different chromatin domain structures within pairwise distance matrices. 

## Requirments:

FISHnet is compatable with Python 3.0.0 or higher. The `requirements.txt` file list all Python libraries that FISHnet depends on, and they can be installed using:

```
pip install -r requirements.txt
```


## Table of Key Parameters:

| Parameter Name  | Type | Description |
| ------------- | ------------- |  ------------- |
| `input_matrix`| array | Pairwise distance matrix from sequential Oligopaints data |
| `distances`  |  list | Ascending distances e.g. [10,20,30,40,...500] that tells FISHnet what thresholds to use. If you want to run FISHnet using only one distance enter a list of one value and keep 'plateau_size' = 0.|
| `plateau_size` | Integer  | The number of adjacent distances required to have the same number of domains for FISHnet to recognize as a stable call. When considering the interplay of `distances` and `plateau_size` keep in mind that having distances that are closely spaced will require a larger plateau size, while having distances that are spaced farther part will require lower a `plateau_size`.|
| `window_size`  | Integer  |  The N by N window used to smooth the pairwise distance matrix. We recommend `window_size` = 2 |
| `merge`  |  Float | Merges Domain calls within N bins of eachother. |
| `size_exclusion`  | Integer | Removes domains less than N bins. |


## Output of FISHnet:

FISHnet outputs two dictionaries that share common keys. The keys correspond to a plateau grouping, where lower key values correspond to inner domain calls, and higher key values are outer domain calls.

| Output Name  | Type | Description |
| ------------- | ------------- |  ------------- |
| `Domains` | dictionary |  A dictionary of domain calls where values correspond to a list of domains stored as tuples. |
| `Distance_scales` | dictionary |  A dictionary of distances where values are a list of distances that were grouped together.|

The output of FISHnet refers to the bins within the `input_matrix` and does not convert them to genomic coordinates.



## Example Code Running FISHnet:

```
from FISHnet_main import FISHnet_main

Domains,Distance_scale  = FISHnet_main(input_matrix=Pairwise_distance,
                          distance= [100,150,200,250,300,350,400,450,500,550,600,650], # distances in nanometers in this example
                          plateau_size=4,
                          window_size=2,
                          size_exclusion=3,
                          merge=3)
```


## Output of Example:

```
Domains

{0: [(0, 4.0),(4.0, 26.0),(26.0, 48.5),(48.5, 65.5),(65.5, 78.5),(78.5, 83)],
 1: [(0, 48.5), (48.5, 65.5), (65.5, 78.5), (78.5, 83)]}
```

In this example two plateau groups were found with the distances used. The domains are stored as a list of tuples indicatig the start and end boundary per each plateau grouping.

```
Distance_scale

{0: [200, 250, 300, 350],
1: [400, 450, 500, 550, 600, 650]}
```
The distances that are grouped together for the plateau groups found in `Domains`. 


## Tutorial Jupyter Notebook

For a jupyter notebook in which you can run FISHnet on HCT116 dataset chr21:34.6-37.1 Mb (Bintu, et. al. Science, 2018). Check out the `FISHnet Tutorial` Folder


 Link to dataset used in the Tutorial: https://github.com/BogdanBintu/ChromatinImaging/tree/master!

