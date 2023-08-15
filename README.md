# NRGA - Non Rigid Gravitational Approach for Point Set Registration
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/static/v1.svg?label=📄%20DOI&message=DOI&color=informational&style=flat-square)](https://ieeexplore.ieee.org/document/8491029)
[![Video](https://img.shields.io/badge/Video-Youtube-ff69b4?style=flat-square)](https://www.youtube.com/watch?v=RrUrjq1uU7o)

This repository presents our work [NRGA](https://www.dfki.de/fileadmin/user_upload/import/9932_NRGA_MainPaper.pdf).

[![Video Thumbnail](media/NRGA_YOUTUBE_IMG.png)](https://www.youtube.com/watch?v=RrUrjq1uU7o) 

## Repository layout

```bash
NRGA/
┣ media/			# media for the README
┣ headers/			# .cfg config files to run NRGA
┣ shaders/                     	# experiment logs are saved here (auto created)
┣ settingsFiles/               	# experiment logs are saved here (auto created)
┣ results/                     	# the checkpoints are saved here (auto created)
┣ libs/                        	# directory to store the data
┃  ┣ ...                       	# detailed instructions in the dataset.md
┣ src/                 	       	# 
┃  ┣ datasets/                 	# class definitions for the datasets
┃  ┗ utils/                    	# helper functions
┗
```

<a name="installation"></a>

### Installation

<a name="conda"></a>

#### Download and Compilation

```bash
git clone --recurse-submodules https://github.com/saali14/NRGA.git 
cd NRGA
mkdir build && cd build
cmake ..
make -j8
```

<span style="color: red">**NOTE**:</span> The code has been tested on Ubuntu 18.04, python 3.8


<a name="cite"></a>

## Cite
If you find this work useful or use any part of the code in this repo, please cite our paper:

```bibtext
@inproceedings{Ali3DV18Nrga,
      title={NRGA: Gravitational Approach for Non-Rigid Point Set Registration}, 
      author={Sk Aziz Ali and Vladislav Golyanik and Didier Stricker},
      booktitle={International Conference on 3D Vision (3DV)},
      year={2018},
      pages={756-765},
      doi={10.1109/3DV.2018.00091}
}
```


