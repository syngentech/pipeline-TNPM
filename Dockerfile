FROM centos
WORKDIR /root
ADD docs /root/docs
Add .tmux.conf /root

# Lock versions
ENV SAMTOOLS_VERSION=1.9 \
    STAR_VERSION=2.6.1a \
    GATK_VERSION=4.0.8.1 \
    GFFREAD_VERSION=0.9.12

# Install softwares
RUN yum update -y && yum install -y \
    bzip2 git gcc java make tmux unzip vim wget zsh \
    bzip2-devel libcurl-devel ncurses-devel openssl-devel xz-devel zlib-devel && \
  sh -c "$(wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O -)"
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
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffread-$GFFREAD_VERSION.Linux_x86_64.tar.gz && \
  tar -xzf gffread-$GFFREAD_VERSION.Linux_x86_64.tar.gz && \
  rm gffread-$GFFREAD_VERSION.Linux_x86_64.tar.gz && \
  ln -s /root/gffread-$GFFREAD_VERSION.Linux_x86_64/gffread /usr/bin/

# Setup environments
WORKDIR /data
CMD ["tmux"]
