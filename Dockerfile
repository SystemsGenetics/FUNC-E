FROM ubuntu:18.04

# install system dependencies
RUN apt-get update -qq && apt-get install -qq -y curl git python3-dev python3-pip && ln -s /usr/bin/python3 python

#changing working directory in Docker container
WORKDIR /app

# copy data from local into Docker container
ADD . /app/

# install python dependencies
RUN pip3 install -r requirements.txt