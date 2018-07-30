FROM broadinstitute/gatk
WORKDIR /

# Lock versions
ENV STAR_VERSION 2.6.1a

# Install softwares
RUN wget https://github.com/alexdobin/STAR/archive/$STAR_VERSION.tar.gz && \
  tar -xzf $STAR_VERSION.tar.gz && \
  rm $STAR_VERSION.tar.gz && \
  mv STAR-$STAR_VERSION star

# Setup environments
ENV PATH /star/bin/Linux_x86_64_static:$PATH
WORKDIR /data
CMD ["bash"]
