# FISHnet - a graph theory approach for calling chromatin domains in sequential Oligopaints data

![FISHnet_video_twitter](https://github.com/user-attachments/assets/e5f2aff0-ad25-4568-83d6-3082c8dc7f60)


## Table of Key Parameters:

| Parameter Name  | Type | Description |
| ------------- | ------------- |  ------------- |
| `Distances`  |  list | Distances i.e. [10,20,30,40,...500] that tells FISHnet what thresholds to use. If you want to run FISHnet using only one distance enter a list of one value and keep 'Plateau Size' = 0.|
| `Plateau Size` | Integer  | The number of adjacent distances required to have the same number of domains for FISHnet to recognize as a stable call. |
| `Smoothing Window`  | Integer  |  The N by N window used to smooth the pairwise distance matrix.  |
| `Merge`  |  Float | Merges Domain calls within N bins of eachother.   |
| `Size Exclusion`  | Integer | Removes domains less than N bins. |

When considering the interplay of 'Distances' and 'Plateau Size' keep in mind that having distances that are closely spaced will require a larger plateau size, while having distances that are spaced farther part will require lower a 'Plateau size'.


## Requirments:

FISHnet runs on Python Version 3 with the following requirments:

contourpy==1.2.1
cycler==0.12.1
fonttools==4.51.0
joblib==1.4.2
kiwisolver==1.4.5
matplotlib==3.8.4
numpy==1.26.4
packaging==24.0
pandas==2.2.2
pillow==10.3.0
pyparsing==3.1.2
python-dateutil==2.9.0.post0
pytz==2024.1
scikit-learn==1.4.2
scipy==1.13.0
six==1.16.0
threadpoolctl==3.5.0
tzdata==2024.1




## Example Code Running FISHnet:




## Output of Example:




## For a jupyter notebook in which you can run FISHnet on HCT116 dataset chr21:34.6-37.1Mb (Bintu, et. al. Science, 2018). Link to dataset: https://github.com/BogdanBintu/ChromatinImaging/tree/master go to tutorials folder! 

