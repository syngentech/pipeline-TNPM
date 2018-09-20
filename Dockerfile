FROM centos
WORKDIR /root
ADD docs /root/docs
ADD pkgs /root/docs
ADD utils /root/utils
Add .tmux.conf /root

# Lock versions
ENV SAMTOOLS_VERSION=1.9 \
    STAR_VERSION=2.6.1a \
    GATK_VERSION=4.0.8.1

# Install softwares
RUN yum update -y && yum install -y \
    bzip2 git gcc java make ruby tcsh tmux unzip vim wget zsh \
    bzip2-devel libcurl-devel ncurses-devel openssl-devel xz-devel zlib-devel && \
  sh -c "$(wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O -)" && \
  curl https://bootstrap.pypa.io/get-pip.py | python && \
  pip install Biopython numpy pandas
RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 && \
  tar -xjf samtools-$SAMTOOLS_VERSION.tar.bz2 && \
  rm samtools-$SAMTOOLS_VERSION.tar.bz2 && \
  cd samtools-$SAMTOOLS_VERSION && \
  ./configure && make && \
  ln -s /root/samtools-$SAMTOOLS_VERSION/samtools /usr/bin/
RUN wget https://github.com/alexdobin/STAR/archive/$STAR_VERSION.tar.gz && \
  tar -xzf $STAR_VERSION.tar.gz && \
  rm $STAR_VERSION.tar.gz && \
  ln -s /root/STAR-$STAR_VERSION/bin/Linux_x86_64_static/STAR /usr/bin/
RUN wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip && \
  unzip gatk-$GATK_VERSION.zip && \
  rm gatk-$GATK_VERSION.zip && \
  ln -s /root/gatk-$GATK_VERSION/gatk /usr/bin/
RUN tar -xzf pkgs/netMHCpan-4.0a.Linux.tar.gz && \
  ln -s /root/netMHCpan-4.0/netMHCpan /usr/bin/ && \
  wget http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.linux.tar.gz && \
  tar -xzf data.linux.tar.gz && \
  rm data.linux.tar.gz && \
  mv data netMHCpan-4.0/Linux_x86_64/

# Setup environments
WORKDIR /data
CMD ["tmux"]
