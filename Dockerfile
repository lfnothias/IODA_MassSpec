FROM jupyter/base-notebook:latest

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gnupg \
    ca-certificates \
    apt-transport-https

# Add any required system dependencies here

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF && \
    echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    mono-devel && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


USER $NB_UID

# Install Python packages from requirements.txt
COPY requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt

COPY pyopenms_wheels /home/jovyan/pyopenms_wheels

# Execute postBuild script
COPY postBuild /tmp/
RUN chmod +x /tmp/postBuild && \
    /tmp/postBuild