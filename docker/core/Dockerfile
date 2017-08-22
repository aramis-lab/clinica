FROM ubuntu:14.04
MAINTAINER GUILLON Jeremy <contact@guillonjeremy.co>

# Build-time variables
# --build-arg
ARG GITLAB_USER
ARG GITLAB_PSWD

RUN apt-get update && apt-get install -y \
    git \
    wget \
    libxml2-dev \
    libxslt-dev \
    gcc

# Miniconda2
RUN wget --progress=bar:force -O /tmp/miniconda2-installer.sh https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda2-installer.sh -b && \
    rm /tmp/miniconda2-installer.sh
ENV PATH="/root/miniconda2/bin:$PATH"

# Clinica dependencies
RUN pip install --upgrade pip
RUN pip install \
    nipype \
    pandas \
    dipy

# Clinica
ENV CLINICA_HOME=/root/clinica
RUN git clone https://$GITLAB_USER:$GITLAB_PSWD@gitlab.icm-institute.org/aramislab/clinica.git $CLINICA_HOME
WORKDIR $CLINICA_HOME
RUN python setup.py install
RUN echo "# Clinica Autocompletion\neval \"\$(register-python-argcomplete clinica)\"\n" >> ~/.bashrc

COPY entrypoint.sh /bin/
RUN chmod +x /bin/entrypoint.sh
# Remove 3 lines in the .bashrc to avoid returning when non-interactive
RUN sed -i '5,7 d' ~/.bashrc
# Force colored prompt: life is better in colors
RUN sed -i 's/#force_color_prompt=yes/force_color_prompt=yes/g' ~/.bashrc
ENTRYPOINT ["/bin/entrypoint.sh"]
