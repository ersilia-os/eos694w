FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit
RUN pip install git+https://github.com/MolecularAI/REINVENT4.git --extra-index-url https://download.pytorch.org/whl/cu113 --extra-index-url https://pypi.anaconda.org/OpenEye/simple

WORKDIR /repo
COPY . /repo
