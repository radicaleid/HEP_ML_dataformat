pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0
pip install -r requirements-cpu.txt  
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.4.0+cpu.html
pip install -e .
python check_acorn.py

lsetup "cmake 3.30.5"
module load gcc/12.3.0