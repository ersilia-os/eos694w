FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia


RUN pip install rdkit-pypi==2022.9.5
RUN pip install protobuf==3.18.3
RUN pip install tensorboardx==2.0
RUN pip install chemprop==1.5.2
RUN pip install tensorboard==2.11.0

RUN pip install git+https://github.com/MolecularAI/REINVENT4@455cd1de86a3b086961497c1aba33488e939b908 --extra-index-url=https://pypi.anaconda.org/OpenEye/simple --extra-index-url https://download.pytorch.org/whl/cu121


WORKDIR /repo
COPY . /repo