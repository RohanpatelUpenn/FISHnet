# FISHnet - a graph theory approach for calling chromatin domains in sequential Oligopaints data

![FISHnet_video_twitter](https://github.com/user-attachments/assets/e5f2aff0-ad25-4568-83d6-3082c8dc7f60)


## Table of Key Parameters:

| Parameter Name  | Type | Discription |
| ------------- | ------------- |  ------------- |
| `Distances`  |  list | Distances i.e. [10,20,30,40,...500] that tells FISHnet what thresholds to use. If you want to run FISHnet using only one distance enter a list of one value and keep 'Plateau Size' = 0.|
| `Plateau Size` | Integer  | The number of adjacent distances required to have the same number of domains for FISHnet to recognize as a stable call. |
| `Smoothing Window`  | Integer  |  The N by N window used to smooth the pairwise distance matrix.  |
| `Merge`  |  Float | Merges Domain calls within N bins of eachother.   |
| `Size Exclusion`  | Integer | Removes domains less than N bins. |

When considering the interplay of 'Distances' and 'Plateau Size' keep in mind that having distances that are closely spaced will require a larger plateau size, while having distances that are spaced farther part will require lower a 'Plateau size'.
