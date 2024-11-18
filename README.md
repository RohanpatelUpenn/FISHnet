# FISHnet - a graph theory approach for calling chromatin domains in sequential Oligopaints data

![FISHnet_video_twitter](https://github.com/user-attachments/assets/e5f2aff0-ad25-4568-83d6-3082c8dc7f60)


## Table of Key Parameters:

| Parameter Name  | Type | Description |
| ------------- | ------------- |  ------------- |
| `input_matrix`| array | Pairwise distance matrix from sequential Oligopaints data |
| `distances`  |  list | Distances i.e. [10,20,30,40,...500] that tells FISHnet what thresholds to use. If you want to run FISHnet using only one distance enter a list of one value and keep 'plateau_size' = 0.|
| `plateau_size` | Integer  | The number of adjacent distances required to have the same number of domains for FISHnet to recognize as a stable call. |
| `window_size`  | Integer  |  The N by N window used to smooth the pairwise distance matrix.  |
| `merge`  |  Float | Merges Domain calls within N bins of eachother.   |
| `size_exclusion`  | Integer | Removes domains less than N bins. |

When considering the interplay of `distances` and `plateau_size` keep in mind that having distances that are closely spaced will require a larger plateau size, while having distances that are spaced farther part will require lower a `plateau_size`.

## Output of FISHnet:

| Output Name  | Type | Description |
| Domains | dictionary |  ------------- |
| Distance_scales | dictionary |  ------------- |


## Requirments:

FISHnet is compatable with Python 3.0.0 or higher. The `requirements.txt` file list all Python libraries that FISHnet depends on, and they can be installed using:

```
pip install -r requirements.txt
```


## Example Code Running FISHnet:




## Output of Example:




## Tutorial Jupyter Notebook

For a jupyter notebook in which you can run FISHnet on HCT116 dataset chr21:34.6-37.1Mb (Bintu, et. al. Science, 2018). Check out the `Tutorial Folder`


 Link to dataset used in the Tutorial: https://github.com/BogdanBintu/ChromatinImaging/tree/master!

