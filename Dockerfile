FROM python:2.7-slim

COPY . /

RUN apt-get update && apt-get install -y \
    apt-utils \
    g++

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

#WORKDIR /
#COPY . /

EXPOSE 7778

CMD ["gunicorn", "--bind", "127.0.0.1:7778", "--workers", "1", "-k", "gevent", "exac:app"]
