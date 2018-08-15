FROM ubuntu:18.04

RUN apt-get update -y --fix-missing
RUN apt-get install git make gcc zlib1g-dev -y

RUN mkdir /hercules
RUN mkdir /input
RUN mkdir /output
WORKDIR /hercules

RUN apt install -y xz-utils
RUN apt install -y g++
RUN git clone https://github.com/BilkentCompGen/hercules.git /hercules
RUN cd src && make

ENV PATH="/hercules/bin:${PATH}"
VOLUME /input
VOLUME /output
ENTRYPOINT ["/hercules/bin/hercules"]
