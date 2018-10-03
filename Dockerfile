FROM sdwfrost/epirecipes:latest

LABEL maintainer="Simon Frost <sdwfrost@gmail.com>"

ENV DEBIAN_FRONTEND noninteractive

USER jovyan

RUN cd $HOME &&
    git clone https://github.com/epirecipes/epicookbook && \
    mv epicookbook/notebooks ${HOME}/notebooks && \
    rm -rf epicookbook

