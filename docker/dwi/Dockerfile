FROM guillonjeremy/clinica:core
MAINTAINER GUILLON Jeremy <contact@guillonjeremy.co>


# ANTs
RUN wget --progress=bar:force -O /tmp/ants.tar.bz2 https://github.com/stnava/ANTs/releases/download/v2.1.0/Linux_Ubuntu14.04.tar.bz2
RUN tar xjvf /tmp/ants.tar.bz2 --strip-components=1 -C /opt/
ENV ANTS_HOME=/opt/ants
ENV PATH=$ANTS_HOME:$PATH


WORKDIR /root
RUN rm -rf /tmp/*
