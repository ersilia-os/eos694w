FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia



RUN pip install rdkit-pypi==2022.9.5
RUN pip install protobuf==3.18.3
RUN pip install tensorboardx==2.0
RUN pip install chemprop==1.5.2
RUN pip install tensorboard==2.11.0

# To install REINVENT4 pip is first clonning it into a temp directory, then running checkout command.
# But Git LFS doesn't able to recognize prior files as pointer. It could be because those prior files
# weren't properly added to Git LFS. Therefore we are disable Git LFS during installation.
ENV GIT_LFS_SKIP_SMUDGE 1

RUN pip install git+https://github.com/MolecularAI/REINVENT4@ee87830e800af22f17b05027bb9378950d5980f7 --extra-index-url=https://pypi.anaconda.org/OpenEye/simple --extra-index-url https://download.pytorch.org/whl/cu113

WORKDIR /repo
COPY . /repo