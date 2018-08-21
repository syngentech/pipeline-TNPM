FROM centos
WORKDIR /root

# Lock versions
ENV STAR_VERSION=2.6.1a \
    GATK_VERSION=4.0.8.1

# Install softwares
RUN yum update -y && \
  yum install -y git java tmux unzip vim wget zsh && \
  sh -c "$(wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O -)"
RUN wget https://github.com/alexdobin/STAR/archive/$STAR_VERSION.tar.gz && \
  tar -xzf $STAR_VERSION.tar.gz && \
  rm $STAR_VERSION.tar.gz && \
  ln -s /root/STAR-2.6.1a/bin/Linux_x86_64_static/STAR /usr/bin/
RUN wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip && \
  unzip gatk-$GATK_VERSION.zip && \
  rm gatk-$GATK_VERSION.zip && \
  ln -s /root/gatk-$GATK_VERSION/gatk /usr/bin/

# Setup environments
WORKDIR /data
ADD README.md /root/
CMD ["tmux"]
