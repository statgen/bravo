FROM python:2.7-slim

COPY . /


RUN apt-get update && apt-get install -y \
    apt-utils \
    autoconf \
    automake \
    bzip2 \
    curl \
    g++ \
    gcc \
    libbz2-dev \
    libcurl4-gnutls-dev \
    liblzma-dev \
    libssl-dev \
    make \
    perl \
    wget \
    zlib1g-dev

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

EXPOSE 80

CMD ["gunicorn", "--bind", "0.0.0.0:80", "--workers", "1", "-k", "gevent", "exac:app"]
