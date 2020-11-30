
FROM gcr.io/oblivion/ubuntu18_py37:latest
 
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install python-minimal python-pip python-tk

RUN python2.7 -m pip install pysam pandas matplotlib seaborn numpy

RUN git clone https://github.com/dayzerodx/counterr.git && cd counterr && python2.7 setup.py install 



