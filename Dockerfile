FROM sdwfrost/epirecipes:latest

LABEL maintainer="Simon Frost <sdwfrost@gmail.com>"

ENV DEBIAN_FRONTEND noninteractive

USER jovyan

RUN cp -r ./notebooks ${HOME}/notebooks

