FROM openms/executables:release-2.5.0

# Installing core dependencies for python
RUN apt-get update -y
RUN apt-get install -y python3-pip python3 build-essential

# Installing Virtual Frame Buffer
RUN apt-get install -y xvfb
RUN apt-get install -y wget
RUN apt-get install -y nodejs

### Mono install
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y gnupg ca-certificates
RUN apt-get -y install sudo 
RUN apt update
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main"
RUN sudo tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt update
RUN apt-get install -y unzip
RUN apt mono-complete
RUN apt update

# Binder Specific Installations
ENV JUPYTER_ENABLE_LAB=yes
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
ENV HOME=/home/${NB_USER}
ENV HOME=/home/jovyan/

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
    
# Installing Requirements
COPY requirements.txt ${HOME}/requirements.txt
RUN pip3 install -r ${HOME}/requirements.txt
RUN pip3 install --no-cache-dir notebook
RUN pip3 install --no-cache-dir jupyterlab
RUN pip3 install --no-cache-dir jupyterlab
RUN pip3 install --no-cache-dir jupyterlab-system-monitor
ENV jupyter lab --notebook-dir=/home/jovyan/

# Copying files
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
