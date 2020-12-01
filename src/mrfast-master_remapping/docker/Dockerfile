FROM ubuntu:16.04
RUN apt-get update -y --fix-missing
RUN apt-get install git make gcc zlib1g-dev -y
RUN mkdir /mrfast
RUN mkdir /input
RUN mkdir /output
WORKDIR /mrfast
RUN git clone --recursive https://github.com/BilkentCompGen/mrfast.git /mrfast
RUN make
ENV PATH="/mrfast:${PATH}"
VOLUME /input
VOLUME /output
ENTRYPOINT ["/mrfast/mrfast"]
