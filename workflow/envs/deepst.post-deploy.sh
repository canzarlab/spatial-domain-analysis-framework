if [ ! -d workflow/scripts/DeepST ] ; then
    git clone https://github.com/JiangBioLab/DeepST workflow/scripts/DeepST
fi

pip3 install torch==1.13.0+cu116 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116  #### GPU
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv torch_geometric -f https://data.pyg.org/whl/torch-1.13.0+cu116.html #### GPU
pip install scikit-network
